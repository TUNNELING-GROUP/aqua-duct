'''
This is driver for aqueduct.
'''

from aqueduct import version as aqueduct_version

from aqueduct.utils import log
from aqueduct.geom.smooth import WindowSmooth

import sys
import os

import ConfigParser

from aqueduct.traj.reader import ReadAmberNetCDFviaMDA

import multiprocessing as mp
import copy
import numpy as np
import cPickle as pickle
import gzip
from aqueduct.traj.paths import GenericPaths, yield_single_paths, InletTypeCodes
from aqueduct.geom.cluster import perform_clustering

from collections import namedtuple

from aqueduct.geom import traces

cpu_count = mp.cpu_count()
optimal_threads = int(1.5 * cpu_count + 1)  # is it really optimal?


__version__ = '20160318 beta'

# optimal_threads = int(2*cpu_count + 1) # is it really optimal?

def CHullCheck(point):
    return CHullCheck.chull.point_within(point)


def CHullCheck_init(args):
    CHullCheck.chull = copy.deepcopy(args[0])


def CHullCheck_pool(chull, threads=optimal_threads):
    return mp.Pool(threads, CHullCheck_init, [(chull,)])


def CHullCheck_exec(chull, points, threads=optimal_threads):
    pool = CHullCheck_pool(chull, threads=threads)
    out = pool.map(CHullCheck, points)
    pool.close()
    pool.join()
    del pool
    return out


class ValveConfig(object):
    def __init__(self):
        self.config = self.get_default_config()

    def __make_options_nt(self, options):
        options_nt = namedtuple('Options', options.keys())
        return options_nt(**options)

    @staticmethod
    def common_config_names():
        # execute - what to do: skip, run
        # load - load previous results form file name
        # save - save results to file name
        return 'execute load save'.split()

    @staticmethod
    def common_traj_data_config_names():
        # scope - scope definition
        # scope_convexhull - take convex hull of scope, true of false
        # object - object definition
        return 'scope scope_convexhull object'.split()

    def global_name(self):
        return 'global settings'

    def cluster_name(self):
        return 'clusterization'
    
    def smooth_name(self):
        return 'smooth'

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
        return self.__make_options_nt(options)

    def get_global_options(self):
        section = self.global_name()
        names = self.config.options(section)
        options = {name: self.config.get(section, name) for name in names}
        return self.__make_options_nt(options)

    def get_stage_options(self, stage):
        assert isinstance(stage, int)
        stage_name = self.stage_names(stage)
        names = self.config.options(stage_name)
        options = {name: self.config.get(stage_name, name) for name in names}
        options.update(self.get_common_traj_data(stage)._asdict())
        return self.__make_options_nt(options)

    def get_cluster_options(self):
        section = self.cluster_name()
        names = self.config.options(section)
        options = {name: self.config.get(section, name) for name in names}
        return self.__make_options_nt(options)

    def get_default_config(self):
        config = ConfigParser.RawConfigParser()

        def common(section):
            for setting in self.common_config_names():
                config.set(section, setting)

        def common_traj_data(section):
            for setting in self.common_traj_data_config_names():
                config.set(section, setting)

        ################
        # global settings
        section = self.global_name()
        config.add_section(section)

        # top - top file name
        # nc - netcdf file name
        config.set(section, 'top')
        config.set(section, 'nc')

        config.set(section, 'pbar', 'tqdm')

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

        config.set(section, 'clear_in_object_info', 'False')

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

        config.set(section, 'apply_smoothing', 'True')

        ################
        # clusterization
        section = self.cluster_name()
        config.add_section(section)
        config.set(section, 'method', 'dbscan')
        config.set(section, 'eps', '5.0')
        config.set(section, 'min_samples', '3')
 
        ################
        # stage V
        # analysis
        section = self.stage_names(4)
        config.add_section(section)

        common(section)

        return config

    def load_config(self, filename):
        self.config.read(filename)

    def save_config_stream(self, fs):
        self.config.write(fs)

    def save_config(self, filename):
        with open(filename, 'w') as fs:
            self.save_config_stream(fs)

def greetings_aqueduct():
    greet = '''------------------------------------------------
           ~ ~ ~ A Q U E D U C T ~ ~ ~
################################################
####        ########        ########        ####
##            ####            ####            ##
#              ##              ##              #
#              ##              ##              #
#              ##              ##              #
#              ##              ##              #
------------------------------------------------'''
    return greet


def sep():
    return '------------------------------------------------'

def asep():
    return '-'*72

def check_res_in_scope(options, scope, res, res_coords):
    if options.scope_convexhull == 'True':
        # find convex hull of protein
        chull = scope.get_convexhull_of_atom_positions()
        current_threads = len(res_coords)
        # current_threads = 1
        if current_threads > optimal_threads:
            current_threads = optimal_threads

        is_res_in_scope = CHullCheck_exec(chull, res_coords, threads=current_threads)
    else:
        res_uids = res.unique_resids()
        res_in_scope_uids = reader.parse_selection(options.scope).unique_resids()
        is_res_in_scope = [ruid in res_in_scope_uids for ruid in res_uids]
    return is_res_in_scope


def get_res_in_scope(is_res_in_scope, res):
    res_new = None
    for iris, r in zip(is_res_in_scope, res.iterate_over_residues()):
        if iris:
            if res_new is None:
                res_new = r
            else:
                res_new += r
    return res_new


if __name__ == "__main__":
    ################################################################################
    # argument parsing
    import argparse
    parser = argparse.ArgumentParser(description="Valve, Aqueduct driver",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--dump-template-config", action="store_true", dest="dump_template_conf", required=False,
                        help="Dumps template config file. Supress all other output or actions.")

    parser.add_argument("-c", action="store", dest="config_file", required=False, help="Config file filename.")
    args = parser.parse_args()
    ################################################################################
    # special option for dumping template config
    config = ValveConfig() # config template
    if args.dump_template_conf:
        import StringIO

        config_dump = StringIO.StringIO()
        config.save_config_stream(config_dump)
        print config_dump.getvalue()
        exit()
    ################################################################################
    # begin!
    log.message(greetings_aqueduct()) # nice greetings
    log.message('Aqueduct version %s' % aqueduct_version())
    log.message('Valve driver version %s' % __version__)
    log.message(sep())
    
    assert args.config_file is not None, "No config file provided."

    # load config file
    with log.fbm('Load configuration file'):
        config.load_config(args.config_file)

    log.message("Optimal threads count: %d" % optimal_threads)

    # get global options
    goptions = config.get_global_options()
    pbar_name = goptions.pbar


    ################################################################################
    # STAGE 0

    with log.fbm('Read trajectory'):
        # read trajectory
        topology = goptions.top
        trajectory = goptions.nc
        reader = ReadAmberNetCDFviaMDA(topology, trajectory)

    ################################################################################
    # STAGE I
    log.message(sep())
    log.message('Starting Stage I: %s' % config.stage_names(0))
    options = config.get_stage_options(0)
    log.message('Execute mode: %s' % options.execute)

    max_frame = reader.number_of_frames - 1
    max_frame = 1000

    # execute?
    if options.execute == 'run':
        # this creates scope

        log.message("Loop over frames - search of residues in object:")
        pbar = log.pbar(max_frame, kind=pbar_name)
        

        scope = reader.parse_selection(options.scope)

        # create some containers
        res_ids_in_object_over_frames = {}
        all_res = None

        for frame in reader.iterate_over_frames():
            if frame > max_frame:
                break
            # current res selection
            res = reader.parse_selection(options.object)

            # check is res are in scope
            if options.scope_convexhull == 'True':
                res_coords = list(res.center_of_mass_of_residues())
                is_res_in_scope = check_res_in_scope(options, scope, res, res_coords)
            else:
                is_res_in_scope = check_res_in_scope(options, scope, res, None)

            # discard res out of scope
            res_new = get_res_in_scope(is_res_in_scope, res)

            # add it to all res in object
            if res_new is not None:
                if all_res:
                    all_res += res_new
                    all_res.uniquify()
                else:
                    all_res = res_new

            # remeber ids of res in object in current frame
            if res_new is not None:
                res_ids_in_object_over_frames.update({frame: res_new.unique_resids()})
            else:
                res_ids_in_object_over_frames.update({frame: []})
            pbar.update(frame)

        pbar.finish()

        log.message("Number of residues to trace: %d" % all_res.unique_resids_number())

        if options.save not in ['None']:
            with log.fbm('Saving data dump in %s file' % options.save):
                with gzip.open(options.save, mode='w', compresslevel=9) as f:
                    pickle.dump({'all_res': all_res, 'res_ids_in_object_over_frames': res_ids_in_object_over_frames}, f)

    elif options.execute in ['skip']:

        if options.load not in ['None']:
            with log.fbm('Loading data dump from %s file' % options.load):
                with gzip.open(options.save, mode='r') as f:
                    loaded_data = pickle.load(f)
                res_ids_in_object_over_frames = {}
                all_res = None
                if 'all_res' in loaded_data.keys():
                    all_res = loaded_data['all_res']
                if 'res_ids_in_object_over_frames' in loaded_data.keys():
                    res_ids_in_object_over_frames = loaded_data['res_ids_in_object_over_frames']

    else:
        raise NotImplementedError('exec mode %s not implemented' % options.execute)

    ################################################################################
    # STAGE II
    log.message(sep())
    log.message('Starting Stage II: %s' % config.stage_names(1))
    options = config.get_stage_options(1)
    log.message('Execute mode: %s' % options.execute)

    # execute?
    if options.execute == 'run':

        if options.clear_in_object_info == 'True':
            log.message('Clear data on residues in object over frames.This will be recalculated.')
            log.message('This will be recalculated on demand.')
            res_ids_in_object_over_frames = {}

        with log.fbm("Init paths container"):
            # type and frames, consecutive elements correspond to residues in all_H2O
            paths = [GenericPaths(resid, min_pf=0, max_pf=max_frame) for resid in all_res.unique_resids()]

        log.message("Trajectory scan:")
        pbar = log.pbar(max_frame, kind=pbar_name)

        scope = reader.parse_selection(options.scope)

        for frame in reader.iterate_over_frames():
            if frame > max_frame:
                break

            all_res_coords = list(all_res.center_of_mass_of_residues())

            # check is res are in scope
            is_res_in_scope = check_res_in_scope(options, scope, all_res, all_res_coords)

            all_resids = [res.first_resid() for res in all_res.iterate_over_residues()]

            for nr, (coord, isscope, resid) in enumerate(zip(all_res_coords, is_res_in_scope, all_resids)):
                if isscope:
                    paths[nr].add_coord(coord)

                    # do we have info on res_ids_in_object_over_frames?
                    if not res_ids_in_object_over_frames.has_key(frame):
                        res = reader.parse_selection(options.object)
                        # discard res out of scope
                        res_new = get_res_in_scope(is_res_in_scope, res)
                        # remeber ids of res in object in current frame
                        if res_new is not None:
                            res_ids_in_object_over_frames.update({frame: res_new.unique_resids()})
                        else:
                            res_ids_in_object_over_frames.update({frame: []})

                    # in scope
                    if resid not in res_ids_in_object_over_frames[frame]:
                        paths[nr].add_scope(frame)
                    else:
                        # in object
                        paths[nr].add_object(frame)

            pbar.update(frame)

        pbar.finish()

        if options.save not in ['None']:
            with log.fbm('Saving data dump in %s file' % options.save):
                with gzip.open(options.save, mode='w', compresslevel=9) as f:
                    pass
                    pickle.dump({'all_res': all_res, 'paths': paths}, f)

    elif options.execute == 'skip':

        if options.load not in ['None']:
            with log.fbm('Loading data dump from %s file' % options.load):
                with gzip.open(options.save, mode='r') as f:
                    loaded_data = pickle.load(f)
                paths = []
                all_res = None
                if 'all_res' in loaded_data.keys():
                    all_res = loaded_data['all_res']
                if 'paths' in loaded_data.keys():
                    paths = loaded_data['paths']

    else:
        raise NotImplementedError('exec mode %s not implemented' % options.execute)


    ################################################################################
    # STAGE III
    log.message(sep())
    log.message('Starting Stage III: %s' % config.stage_names(2))
    options = config.get_stage_options(2)
    log.message('Execute mode: %s' % options.execute)

    # execute?
    if options.execute == 'run':

        if options.discard_empty_paths == 'True':
            with log.fbm("Discard residues with empty paths"):
                # zipped = [(wat,path) for wat,path in zip(all_H2O.iterate_over_residues(),paths) if len(path.frames) > 0]
                all_res_, paths_ = zip(
                    *[(r, path) for r, path in zip(all_res.iterate_over_residues(), paths) if len(path.frames) > 0])
                all_res = all_res_
                paths = paths_
                del all_res_
                del paths_

        log.message("Create separate paths:")

        pbar = log.pbar(len(paths), kind=pbar_name)
        spaths = [sp for sp, nr in yield_single_paths(paths, progress=True) if pbar.update(nr + 1) is None]


        pbar.finish()

        if options.save not in ['None']:
            with log.fbm('Saving data dump in %s file' % options.save):
                with gzip.open(options.save, mode='w', compresslevel=9) as f:
                    pickle.dump({'paths': paths, 'spaths': spaths}, f)

    elif options.execute == 'skip':

        if options.load not in ['None']:
            with log.fbm('Loading data dump from %s file' % options.load):
                with gzip.open(options.save, mode='r') as f:
                    loaded_data = pickle.load(f)
                spaths = []
                if 'spaths' in loaded_data.keys():
                    spaths = loaded_data['spaths']
                if 'paths' in loaded_data.keys():
                    paths = loaded_data['paths']

    else:
        raise NotImplementedError('exec mode %s not implemented' % options.execute)

    ################################################################################
    # STAGE IV
    log.message(sep())
    log.message('Starting Stage IV: %s' % config.stage_names(3))
    options = config.get_stage_options(3)
    log.message('Execute mode: %s' % options.execute)
    coptions = config.get_cluster_options()
    
    # execute?
    if options.execute == 'run':
        # apply smoothing?
        if options.apply_smoothing == 'True':
            with log.fbm('Applying smoothing'):
                smooth = WindowSmooth() # TODO: read smooth options
                for sp in spaths:
                    sp.apply_smoothing(smooth)
         
        # find coords of inlets, type and id
        inlet_coords = []
        inlet_type = []
        inlet_id = []
        inlet_spnr = []
        for nr,sp in enumerate(spaths):
            for inlet,itype in sp.coords_filo:
                inlet_coords.append(inlet.tolist())
                inlet_type.append(itype)
                inlet_id.append(sp.id)
                inlet_spnr.append(nr)
                #print nr,sp.id,inlet,itype
        
        # find nr of inlets
        nr_of_inlets = len(inlet_id)

        if nr_of_inlets > 0:
            with log.fbm("Performing clusterization"):
                assert coptions.method == 'dbscan', 'Unknown clusterization method %s.' % coptions.method
                # perform clusterization
                clusters = perform_clustering(inlet_coords, eps=float(coptions.eps), min_samples=int(coptions.min_samples))
        else:
            log.message("No inlets found. Clusterization skipped.")
            clusters = np.array([])

        if options.save not in ['None']:
            with log.fbm('Saving data dump in %s file' % options.save):
                with gzip.open(options.save, mode='w', compresslevel=9) as f:
                    pickle.dump({'inlet_coords':inlet_coords,
                                 'inlet_type':inlet_type,
                                 'inlet_id':inlet_id,
                                 'inlet_spnr':inlet_spnr,
                                 'clusters':clusters},f)

    elif options.execute == 'skip':

        if options.load not in ['None']:
            with log.fbm('Loading data dump from %s file' % options.load):
                with gzip.open(options.save, mode='r') as f:
                    loaded_data = pickle.load(f)
                    inlet_coords = []
                    inlet_type = []
                    inlet_id = []
                    inlet_spnr = []
                    clusters = np.array([])
                    if 'inlet_coords' in loaded_data.keys():
                        inlet_coords = loaded_data['inlet_coords']
                    if 'inlet_type' in loaded_data.keys():
                        coords_inlets = loaded_data['inlet_type']
                    if 'inlet_id' in loaded_data.keys():
                        inlet_id = loaded_data['inlet_id']
                    if 'inlet_spnr' in loaded_data.keys():
                        inlet_spnr = loaded_data['inlet_spnr']
                    if 'clusters' in loaded_data.keys():
                        clusters = loaded_data['clusters']
    else:
        raise NotImplementedError('exec mode %s not implemented' % options.execute)

    ################################################################################
    # STAGE V
    log.message(sep())
    log.message('Starting Stage V: %s' % config.stage_names(4))
    options = config.get_stage_options(4)
    log.message('Execute mode: %s' % options.execute)

    # execute?
    if options.execute == 'run':
        # file hSobrietyandle?
        if options.save not in ['None']:
            fh = open(options.save,'w')
            log.message('Using user provided file.')
            log.message(sep())
            log.message('')
        else:
            log.message('Using standard output.')
            log.message(sep())
            log.message('')
            fh = sys.stdout

        ############
        print >> fh, asep()
        print >> fh, 'Aqueduct analysis'

        ############
        print >> fh, asep()
        print >> fh, "Number of traceable residues:",len(paths)
        print >> fh, "Number of separate paths:",len(spaths)

        ############
        print >> fh, asep()
        print >> fh, "List of separate paths"
        header_template = " ".join(['%7s']*7)
        print >> fh, header_template % tuple("Nr ID Begin INP OBJ OUT End".split())
        for nr,sp in enumerate(spaths):
            line = []
            for e in (nr,sp.id,sp.begins,len(sp.path_in),len(sp.path_object),len(sp.path_out), sp.ends):
                line += ["%7d" % e]
            print >> fh, " ".join(line)

        ############
        print >> fh, asep()
        print >> fh, "Spearate paths lenghts"
        header_template = " ".join(['%7s']*2+['%9s']*3)
        print >> fh, header_template % tuple("Nr ID INP OBJ OUT".split())
        line_template = " ".join(['%7d']*2 + ['%9.1f']*3)
        for nr,sp in enumerate(spaths):
            line = [nr, sp.id]
            for e in sp.coords_in,sp.coords_object,sp.coords_out:
                if len(e) > 1:
                    line += [sum(traces.diff(e))]
                else:
                    line += [float('nan')]
            print >> fh, line_template % tuple(line)

        ############
        print >> fh, asep()
        print >> fh, "Spearate paths average step lenghts"
        header_template = " ".join(['%7s']*2+['%8s']*6)
        print >> fh, header_template % tuple("Nr ID INP INPstd OBJ OBJstd OUT OUTstd".split())
        line_template = " ".join(['%7d']*2 + ['%8.2f','%8.3f']*3)
        for nr,sp in enumerate(spaths):
            line = [nr, sp.id]
            for e in sp.coords_in,sp.coords_object,sp.coords_out:
                if len(e) > 1:
                    line += [np.mean(traces.diff(e))]
                    line += [np.std(traces.diff(e))]
                else:
                    line += [float('nan')]
                    line += [float('nan')]
            print >> fh, line_template % tuple(line)

        ############
        print >> fh, asep()
        print >> fh, "Number of inlets:",len(inlet_id)
        no_of_clusters = len(set([c for c in clusters if c != -1]))
        print >> fh, "Number of clusters:", no_of_clusters
        print >> fh, "Outliers:",{True:'yes',False:'no'}[-1 in clusters]
        
        ############
        print >> fh, asep()
        print >> fh, "Clusters summary"
        clusters_list = list(set(clusters.tolist()))
        clusters_list.sort()
        header_template = " ".join(['%7s']*3)
        print >> fh, header_template % tuple("Nr Cluster Size".split())
        line_template = " ".join(['%7d']*3)
        for nr,c in enumerate(clusters_list):
            line = [nr,int(c),clusters.tolist().count(c)]
            print >> fh, line_template % tuple(line)
        
        ############
        print >> fh, asep()
        print >> fh, "Spearate paths inlets clusters"
        header_template = " ".join(['%7s']*2+['%7s']*2+['%8s']*1)
        print >> fh, header_template % tuple("Nr ID INP OUT Type".split())
        line_template = " ".join(['%7d']*2 + ['%7s']*2+ ['%8s']*1)
        for nr,sp in enumerate(spaths):
            line = [nr, sp.id]
            ids = [nr for nr,iid in enumerate(inlet_id) if iid == sp.id]
            inp,out = None,None
            if len(ids) > 0:
                for iid in ids:
                    if inlet_type[iid] == InletTypeCodes.inlet_in_code:
                        inp = clusters[iid]
                    if inlet_type[iid] == InletTypeCodes.inlet_out_code:
                        out = clusters[iid]
            line.extend([str(inp),str(out)])
            # type?
            if inp is None and out is None:
                ctype = 'internal'
            if inp is not None and out is None:
                ctype = 'incoming'
            if inp is None and out is not None:
                ctype = 'outgoing'
            if inp is not None and out is not None:
                if inp == out:
                    ctype = 'in-out'
                else:
                    ctype = 'through'
            line.append(ctype)
            print >> fh, line_template % tuple(line)

        if options.save not in ['None']:
            fh.close()
    elif options.execute == 'skip':
        pass
    else:
        raise NotImplementedError('exec mode %s not implemented' % options.execute)
