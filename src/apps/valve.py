'''
This is driver for aqueduct.
'''

import ConfigParser

from aqueduct.traj.reader import ReadAmberNetCDFviaMDA
from aqueduct.utils import log

import multiprocessing as mp
import copy
import numpy as np

cpu_count = mp.cpu_count()
optimal_threads = int(1.5*cpu_count + 1) # is it really optimal?
#optimal_threads = int(2*cpu_count + 1) # is it really optimal?

def CHullCheck(point):
    return CHullCheck.chull.point_within(point)

def CHullCheck_init(args):
    CHullCheck.chull = copy.deepcopy(args[0])

def CHullCheck_pool(chull,threads=optimal_threads):
    return mp.Pool(threads, CHullCheck_init, [(chull,)])

def CHullCheck_exec(chull,points,threads=optimal_threads):
    pool = CHullCheck_pool(chull,threads=threads)
    out = pool.map(CHullCheck, points)
    pool.close()
    pool.join()
    del pool
    return out


class ValveConfig(object):
    def __init__(self):
        self.config = self.get_default_config()

    @staticmethod
    def common_config_names():
        # exec - what to do: skip, run
        # load - load previous results form file name
        # save - save results to file name
        return 'exec load save'.split()

    @staticmethod
    def common_traj_data_config_names():
        # top - top file name
        # nc - netcdf file name
        # scope - scope definition
        # scope_convexhull - take convex hull of scope, true of false
        # object - object definition
        return 'top nc scope scope_convexhull object'.split()

    def stage_names(self, nr=None):
        if nr is None:
            return [self.stage_names(nr) for nr in range(5)]
        else:
            if nr == 0:
                return 'find traceable residues'
            elif nr == 1:
                return 'find raw paths'
            elif nr == 2:
                return 'create separate frames'
            elif nr == 3:
                return 'inlets clusterisation'
            elif nr == 4:
                return 'analysis'
        raise NotImplementedError('Stage %r is not implemented.' % nr)

    def get_common_traj_data(self, stage):
        assert isinstance(stage, int)
        options = {name: None for name in self.common_traj_data_config_names()}
        for nr in range(stage + 1)[::-1]:
            section = self.stage_names(nr)
            for name in self.common_traj_data_config_names():
                if self.config.has_option(section, name):
                    value = self.config.get(section, name)
                    if (value is not None) and (options[name] is None):
                        options.update({name: value})
        return options

    def get_options(self, stage):
        assert isinstance(stage, int)
        stage_name = self.stage_names(stage)
        names = self.config.options(stage_name)
        options = {name: self.config.get(stage_name, name) for name in names}
        options.update(self.get_common_traj_data(stage))
        return options

    def get_default_config(self):
        config = ConfigParser.RawConfigParser()

        def common(section):
            for setting in self.common_config_names():
                config.set(section, setting)

        def common_traj_data(section):
            for setting in self.common_traj_data_config_names():
                config.set(section, setting)

        ################
        # stage I
        # find traceable residues
        section = self.stage_names(0)
        config.add_section(section)

        common(section)
        common_traj_data(section)

        ################
        # stage II
        # find raw paths
        section = self.stage_names(1)
        config.add_section(section)

        common(section)
        common_traj_data(section)


        ################
        # stage III
        # create separate frames
        section = self.stage_names(2)
        config.add_section(section)

        common(section)

        config.set(section, 'discard_empty_paths', 'True')

        ################
        # stage IV
        # inlets clusterisation
        section = self.stage_names(3)
        config.add_section(section)

        common(section)

        ################
        # stage V
        # analysis
        section = self.stage_names(4)
        config.add_section(section)

        common(section)

        return config

    def load_config(self, filename):
        self.config.read(filename)

    def save_config_stream(self,fs):
        self.config.write(fs)

    def save_config(self, filename):
        with open(filename, 'w') as fs:
            self.save_config_stream(fs)


if __name__ == "__main__":
    # argument parsing
    import argparse

    parser = argparse.ArgumentParser(description="Valve, Aqueduct driver",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--dump-template-config", action="store_true", dest="dump_template_conf", required=False, help="Dumps template config file. Supress all other output or actions.")

    parser.add_argument("-c", action="store", dest="config_file", required=False, help="Config file filename.")

    args = parser.parse_args()

    config = ValveConfig()

    # special option for dumping template config

    if args.dump_template_conf:
        import StringIO
        config_dump = StringIO.StringIO()
        config.save_config_stream(config_dump)
        print config_dump.getvalue()
        exit()

    assert args.config_file is not None

    # load config file
    log.message('Load config file')
    config.load_config(args.config_file)

    # prepare traj_dat_options
    top = None
    nc = None
    traj_scope = None
    #traj_object = None

    scope_convexhull = None

    reader = None

    # STAGE I
    log.message('Starting Stage I')
    options = config.get_options(0)
    print options
    # execute?
    if options['exec'] == 'run':
        log.message('exec mode is run')

        # this reads trajectory
        log.message('Try to read trajectory')
        if not (top == options['top'] and nc == options['nc']):
            top,nc = options['top'],options['nc']
            reader = ReadAmberNetCDFviaMDA(top, nc)
            max_frame = reader.number_of_frames - 1
        # this creates scope
        log.message('Create scope')
        if traj_scope != options['scope']:
            traj_scope = options['scope']
            scope = reader.parse_selection(traj_scope)

        traj_object = options['object']
        scope_convexhull = options['scope_convexhull']

        log.message("Loop over frames: search of waters in object...")
        pbar = log.pbar(max_frame)

        # create some containers
        res_ids_in_object_over_frames = {}
        all_res = None

        max_frame = 100
        for frame in reader.iterate_over_frames():
            if frame > max_frame:
                break
            # current res selection
            res = reader.parse_selection(traj_object)


            # to check if res is in scope we have two methods:
            # is it within convex hull of scope
            # is it in scope

            if scope_convexhull == 'True':
                res_coords = list(res.center_of_mass_of_residues())
                # find convex hull of protein
                chull = scope.get_convexhull_of_atom_positions()
                current_threads = len(res_coords)
                #current_threads = 1
                if current_threads > optimal_threads:
                    current_threads = optimal_threads

                is_res_in_scope = CHullCheck_exec(chull, res_coords, threads=current_threads)
            else:
                res_uids = res.unique_resids()
                res_in_scope_uids = reader.parse_selection(traj_scope).unique_resids()
                is_res_in_scope = [ruid in res_in_scope_uids for ruid in res_uids]

            # discard res out of scope
            res_new = None
            for iris,r in zip(is_res_in_scope,res.iterate_over_residues()):
                if iris:
                    if res_new is None:
                        res_new = r
                    else:
                        res_new += r


            # add it to all res in object
            if res_new is not None:
                if all_res:
                    all_res += res_new
                    all_res.uniquify()
                else:
                    all_res = res_new

            print all_res

            # remeber ids of res in object in current frame
            if res_new is not None:
                res_ids_in_object_over_frames.update({frame:res_new.unique_resids()})
            else:
                res_ids_in_object_over_frames.update({frame:[]})
            pbar.update(frame)

        pbar.finish()

        log.message("Number of residues to trace: %d" % all_res.unique_resids_number())

    elif options['exec'] == 'skip':
        log.message('exec mode is skip')
    else:
        raise NotImplementedError('exec mode %s not implemented' % options['exec'])
















