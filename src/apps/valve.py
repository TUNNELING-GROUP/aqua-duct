#!/bin/env python2

"""
This is driver for aqueduct.
"""

import ConfigParser
import cPickle as pickle
import copy
import gzip
import multiprocessing as mp
import os
import sys
from collections import namedtuple, OrderedDict
import shlex
import roman

import MDAnalysis as mda
import numpy as np

from aqueduct import greetings as greetings_aqueduct
from aqueduct import version as aqueduct_version
from aqueduct import version_nice as aqueduct_version_nice
from aqueduct.geom import traces
from aqueduct.geom.cluster import perform_clustering
from aqueduct.geom.smooth import WindowSmooth, MaxStepSmooth, WindowOverMaxStepSmooth, ActiveWindowSmooth, \
    ActiveWindowOverMaxStepSmooth, DistanceWindowSmooth, DistanceWindowOverMaxStepSmooth
from aqueduct.traj.dumps import TmpDumpWriterOfMDA
from aqueduct.traj.paths import GenericPaths, yield_single_paths, InletTypeCodes
from aqueduct.traj.reader import ReadAmberNetCDFviaMDA
from aqueduct.traj.selections import CompactSelectionMDA, SelectionMDA
from aqueduct.utils import log
from aqueduct.utils.helpers import range2int

# TODO: Move it to separate module
cpu_count = mp.cpu_count()
# optimal_threads = int(2*cpu_count + 1) # is it really optimal?
optimal_threads = int(1.5 * cpu_count + 1)  # is it really optimal?


def version():
    return 0, 4, 1


def version_nice():
    return '.'.join(map(str, version()))


__mail__ = 'Tomasz Magdziarz <tomasz.magdziarz@polsl.pl>'
__version__ = version_nice()


###############################################################################
# configuration file helpers

class ConfigSpecialNames:
    special_names_dict = {'none': None,
                          'null': None,
                          'true': True,
                          'false': False,
                          'yes': True,
                          'no': False}

    def special_name(self, name):
        if isinstance(name, (str, unicode)):
            if name.lower() in self.special_names_dict:
                return self.special_names_dict[name.lower()]
        return name


class ValveConfig(object, ConfigSpecialNames):
    def __init__(self):
        self.config = self.get_default_config()

    def __make_options_nt(self, input_options):
        # options = {opt: self.special_name(input_options[opt]) for opt in input_options}
        options = dict(
            ((opt, self.special_name(input_options[opt])) for opt in input_options))  # This is due to old pydev
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

    @staticmethod
    def global_name():
        return 'global'

    @staticmethod
    def cluster_name():
        return 'clusterization'

    @staticmethod
    def smooth_name():
        return 'smooth'

    def stage_names(self, nr=None):
        if nr is None:
            return [self.stage_names(nr) for nr in range(5)]
        else:
            if nr == 0:
                return 'traceable_residues'
            elif nr == 1:
                return 'raw_paths'
            elif nr == 2:
                return 'separate_paths'
            elif nr == 3:
                return 'inlets_clusterization'
            elif nr == 4:
                return 'analysis'
            elif nr == 5:
                return 'visualize'
        raise NotImplementedError('Stage %r is not implemented.' % nr)

    def get_common_traj_data(self, stage):
        assert isinstance(stage, int)
        # options = {name: None for name in self.common_traj_data_config_names()}
        options = dict(((name, None) for name in self.common_traj_data_config_names()))
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
        # options = {name: self.config.get(section, name) for name in names}
        options = dict(((name, self.config.get(section, name)) for name in names))
        return self.__make_options_nt(options)

    def get_stage_options(self, stage):
        assert isinstance(stage, int)
        stage_name = self.stage_names(stage)
        names = self.config.options(stage_name)
        # options = {name: self.config.get(stage_name, name) for name in names}
        options = dict(((name, self.config.get(stage_name, name)) for name in names))
        if stage in [0, 1]:
            options.update(self.get_common_traj_data(stage)._asdict())
        return self.__make_options_nt(options)

    def get_cluster_options(self):
        section = self.cluster_name()
        names = self.config.options(section)
        # options = {name: self.config.get(section, name) for name in names}
        options = dict(((name, self.config.get(section, name)) for name in names))
        return self.__make_options_nt(options)

    def get_smooth_options(self):
        section = self.smooth_name()
        names = self.config.options(section)
        # options = {name: self.config.get(section, name) for name in names}
        options = dict(((name, self.config.get(section, name)) for name in names))
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

        config.set(section, 'pbar', 'simple')

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
        config.set(section, 'sort_by_id', 'True')
        config.set(section, 'apply_smoothing', 'False')
        config.set(section, 'apply_soft_smoothing', 'True')

        ################
        # stage IV
        # inlets clusterisation
        section = self.stage_names(3)
        config.add_section(section)

        common(section)

        ################
        # smooth
        section = self.smooth_name()
        config.add_section(section)
        config.set(section, 'method', 'window')
        # config.set(section, 'recursive', '0')
        # config.set(section, 'window', '5')
        # config.set(section, 'function', 'mean')

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
        config.remove_option(section, 'load')

        config.set(section, 'dump_config', 'True')

        ################
        # stage VI
        # visualize
        section = self.stage_names(5)
        config.add_section(section)
        common(section)
        config.remove_option(section, 'load')
        config.remove_option(section, 'save')

        # visualize spaths, all paths in one object
        config.set(section, 'all_paths_raw', 'True')
        config.set(section, 'all_paths_smooth', 'True')
        config.set(section, 'all_paths_split', 'True')  # split by in obj out
        config.set(section, 'all_paths_raw_io', 'True')
        config.set(section, 'all_paths_smooth_io', 'True')

        # visualize spaths, separate objects
        config.set(section, 'paths_raw', 'True')
        config.set(section, 'paths_smooth', 'True')
        config.set(section, 'paths_states', 'True')
        config.set(section, 'paths_raw_io', 'True')
        config.set(section, 'paths_smooth_io', 'True')

        # visualize clusters
        config.set(section, 'inlets_clusters', 'True')

        # show protein
        config.set(section, 'show_molecule', 'None')
        config.set(section, 'show_molecule_frames', '0')

        return config

    def load_config(self, filename):
        self.config.read(filename)

    def save_config_stream(self, fs):
        self.config.write(fs)

    def save_config(self, filename):
        with open(filename, 'w') as fs:
            self.save_config_stream(fs)

    def dump_config(self):
        output = []
        options = [self.get_global_options()] + \
                  [self.get_stage_options(stage) for stage in range(5)] + \
                  [self.get_cluster_options(), self.get_smooth_options()]
        names = [self.global_name()] + \
                [self.stage_names(stage) for stage in range(5)] + \
                [self.cluster_name(), self.smooth_name()]

        for o, n in zip(options, names):
            output.append('[%s]' % n)
            for k in o._asdict().keys():
                output.append('%s = %s' % (k, str(o._asdict()[k])))

        return output


################################################################################
# convex hull helpers
# TODO: Move it to separate module

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


################################################################################
# in scope helpers

def check_res_in_scope(options, scope, res, res_coords):
    if options.scope_convexhull:
        if len(res_coords) == 0:
            return []
        # find convex hull of protein
        chull = scope.get_convexhull_of_atom_positions()
        current_threads = len(res_coords)
        # current_threads = 1
        if current_threads > optimal_threads:
            current_threads = optimal_threads

        is_res_in_scope = CHullCheck_exec(chull, res_coords, threads=current_threads)
    else:
        if res.unique_resids_number() == 0:
            return []
        res_in_scope_uids = reader.parse_selection(options.scope).unique_resids(ikwid=True)
        is_res_in_scope = [r.unique_resids(ikwid=True) in res_in_scope_uids for r in res.iterate_over_residues()]
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


################################################################################
# separators - logging

def sep():
    return log.gsep(sep='-', times=48)


def asep():
    return log.gsep(sep='=', times=72)


################################################################################
# save - load helpers

def save_dump(filename, data_to_save, **kwargs):
    with log.fbm('Saving data dump in %s file' % filename):
        with gzip.open(filename, mode='w', compresslevel=9) as f:
            # first version:
            pickle.dump({'version': version(),
                         'aqueduct_version': aqueduct_version()}, f)
            # then data to save:
            pickle.dump(data_to_save, f)
            # then other kwargs
            pickle.dump(kwargs, f)


def load_dump(filename):
    with log.fbm('Loading data dump from %s file' % filename):
        with gzip.open(filename, mode='r') as f:
            # version!
            loaded_data = pickle.load(f)
            check_versions(loaded_data)
            # loaded data!
            loaded_data = pickle.load(f)
            # enything else?
            try:
                other_data = pickle.load(f)
                loaded_data.update({'other_data': other_data})
            except:
                pass
        return loaded_data
        # loaded_data_nt = namedtuple('LoadedData', loaded_data.keys())
        # return loaded_data_nt(**loaded_data)


def check_version_compilance(current, loaded, what):
    if current[0] > loaded[0]:
        log.error(
            'Loaded data has %s major version lower then the application, possible problems with API compilance.' % what)
    if current[0] < loaded[0]:
        log.error(
            'Loaded data has %s major version higher then the application, possible problems with API compilance.' % what)
    if current[1] > loaded[1]:
        log.warning(
            'Loaded data has %s minor version lower then the application, possible problems with API compilance.' % what)
    if current[1] < loaded[1]:
        log.warning(
            'Loaded data has %s minor version higher then the application, possible problems with API compilance.' % what)


def check_versions(version_dict):
    assert isinstance(version_dict, (dict, OrderedDict)), "File is corrupted, cannot read version data."
    assert 'version' in version_dict, "File is corrupted, cannot read version data."
    assert 'aqueduct_version' in version_dict, "File is corrupted, cannot read version data."
    check_version_compilance(aqueduct_version(), version_dict['aqueduct_version'], 'Aqueduct')
    check_version_compilance(version(), version_dict['version'], 'Valve')


################################################################################
# save - load per stage

def save_stage_dump(name, **kwargs):
    # check if name is None
    if name is not None:
        # ok, check if some of kwargs have to be changed
        data_to_save = {}
        for key, value in kwargs.iteritems():
            # SelectionMDA
            if isinstance(value, SelectionMDA):
                value = CompactSelectionMDA(value)
            # options
            if 'options' in key:
                value = value._asdict()
            data_to_save.update({key: value})
            # now we are redy to save
        save_dump(name, data_to_save)


def load_stage_dump(name, reader=None):
    if name is not None:
        loaded_data = {}
        for key, value in load_dump(name).iteritems():
            # CompactSelectionMDA
            if isinstance(value, CompactSelectionMDA):
                value = value.toSelectionMDA(reader)
            loaded_data.update({key: value})
        return loaded_data


################################################################################

def get_smooth_method(soptions):
    assert soptions.method in ['window', 'mss', 'window_mss', 'awin',
                               'awin_mss', 'dwin', 'dwin_mss'], 'Unknown smoothing method %s.' % soptions.method

    opts = {}
    if 'recursive' in soptions._asdict():
        opts.update({'recursive': int(soptions.recursive)})

    def window_opts():
        if 'window' in soptions._asdict():
            opts.update({'window': int(float(soptions.window))})
        function_opts()

    def awin_dwin_opts():
        if 'window' in soptions._asdict():
            opts.update({'window': float(soptions.window)})
        function_opts()

    def function_opts():
        if 'function' in soptions._asdict():
            assert soptions.function in ['mean', 'median'], 'Unknown smoothing function %s.' % soptions.function
            if soptions.function == 'mean':
                opts.update({'function': np.mean})
            if soptions.function == 'median':
                opts.update({'function': np.median})

    def mss_opts():
        if 'step' in soptions._asdict():
            opts.update({'step': float(soptions.step)})

    if soptions.method == 'window':
        window_opts()
        smooth = WindowSmooth(**opts)
    elif soptions.method == 'awin':
        awin_dwin_opts()
        smooth = ActiveWindowSmooth(**opts)
    elif soptions.method == 'dwin':
        awin_dwin_opts()
        smooth = DistanceWindowSmooth(**opts)
    elif soptions.method == 'mss':
        mss_opts()
        smooth = MaxStepSmooth(**opts)
    elif soptions.method == 'window_mss':
        window_opts()
        mss_opts()
        smooth = WindowOverMaxStepSmooth(**opts)
    elif soptions.method == 'awin_mss':
        awin_dwin_opts()
        mss_opts()
        smooth = ActiveWindowOverMaxStepSmooth(**opts)
    elif soptions.method == 'dwin_mss':
        awin_dwin_opts()
        mss_opts()
        smooth = DistanceWindowOverMaxStepSmooth(**opts)

    return smooth


################################################################################

def valve_begin():
    log.message(greetings_aqueduct())  # nice greetings
    log.message('Aqueduct version %s' % aqueduct_version_nice())
    log.message('Valve driver version %s' % version_nice())
    log.message(sep())


def valve_end():
    log.message(sep())
    log.message('Let the Valve be always open!')
    log.message('Goodby!')


def valve_load_config(filename, config):
    assert filename is not None, "No config file provided."
    with log.fbm('Load configuration file'):
        config.load_config(filename)


def valve_read_trajectory(top, traj):
    with log.fbm('Read trajectory'):
        # read trajectory
        traj_list = shlex.split(traj)
        return ReadAmberNetCDFviaMDA(top, traj_list)
        # reader = ReadDCDviaMDA(topology, trajectory)


def valve_begin_stage(stage, config):
    log.message(sep())
    log.message('Starting Stage %s: %s' % (roman.toRoman(stage + 1), config.stage_names(stage)))
    options = config.get_stage_options(stage)
    log.message('Execute mode: %s' % options.execute)
    return options


def valve_exec_stage(stage, config, stage_run, reader=None, no_io=False,
                     **kwargs):
    options = valve_begin_stage(stage, config)

    # execute?
    if options.execute in ['run']:
        result = stage_run(config, options, reader=reader, **kwargs)
        if not no_io:
            ###########
            # S A V E #
            ###########
            save_stage_dump(options.save, **result)
    elif options.execute in ['skip']:
        if not no_io:
            ###########
            # L O A D #
            ###########
            if options.load:
                result = load_stage_dump(options.load, reader=reader)
    else:
        raise NotImplementedError('exec mode %s not implemented' % options.execute)
    # remove options stuff
    if result is not None:
        return dict(((key, val) for key, val in result.iteritems() if 'options' not in key))


################################################################################
# stages run

def stage_I_run(config, options,
                reader=None,
                max_frame=None,
                **kwargs):
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
        if options.scope_convexhull:
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
            res_ids_in_object_over_frames.update({frame: res_new.unique_resids(ikwid=True)})
        else:
            res_ids_in_object_over_frames.update({frame: []})
        pbar.update(frame)

    pbar.finish()

    if all_res is None:
        raise ValueError("No traceable residues was found.")

    log.message("Number of residues to trace: %d" % all_res.unique_resids_number())

    return {'all_res': all_res,
            'res_ids_in_object_over_frames': res_ids_in_object_over_frames,
            'options': options}


################################################################################

def stage_II_run(config, options,
                 reader=None,
                 all_res=None,
                 res_ids_in_object_over_frames=None,
                 max_frame=None,
                 **kwargs):
    if options.clear_in_object_info:
        log.message('Clear data on residues in object over frames.')
        log.message('This will be recalculated on demand.')
        res_ids_in_object_over_frames = {}

    with log.fbm("Init paths container"):
        paths = dict(((resid, GenericPaths(resid, min_pf=0, max_pf=max_frame)) for resid in
                      all_res.unique_resids(ikwid=True)))

    log.message("Trajectory scan:")
    pbar = log.pbar(max_frame, kind=pbar_name)

    scope = reader.parse_selection(options.scope)

    for frame in reader.iterate_over_frames():
        if frame > max_frame:
            break

        all_res_coords = list(all_res.center_of_mass_of_residues())  # this uses iterate over residues

        # check if is res are in scope
        is_res_in_scope = check_res_in_scope(options, scope, all_res, all_res_coords)

        all_resids = [res.first_resid() for res in all_res.iterate_over_residues()]

        for nr, (coord, isscope, resid) in enumerate(zip(all_res_coords, is_res_in_scope, all_resids)):
            # the point is that nr is not pointing to correct element in paths
            # now, nr is useless because paths is a dictionary, use resids instead
            assert paths[resid].id == resid, \
                "Internal error. Paths IDs not synced with resids. \
                 Please send a bug report to developer(s): %s" % __mail__
            if isscope:
                paths[resid].add_coord(coord)

                # do we have info on res_ids_in_object_over_frames?
                if frame not in res_ids_in_object_over_frames:
                    res = reader.parse_selection(options.object)
                    # discard res out of scope
                    res_new = get_res_in_scope(is_res_in_scope, res)
                    # remeber ids of res in object in current frame
                    if res_new is not None:
                        res_ids_in_object_over_frames.update({frame: res_new.unique_resids(ikwid=True)})
                    else:
                        res_ids_in_object_over_frames.update({frame: []})

                # in scope
                if resid not in res_ids_in_object_over_frames[frame]:
                    paths[resid].add_scope(frame)
                else:
                    # in object
                    paths[resid].add_object(frame)

        pbar.update(frame)

    pbar.finish()

    return {'all_res': all_res, 'paths': paths, 'options': options}


################################################################################

def stage_III_run(config, options,
                  paths=None,
                  **kwargs):
    soptions = config.get_smooth_options()

    if options.discard_empty_paths:
        with log.fbm("Discard residues with empty paths"):
            for key in paths.keys():
                if len(paths[key].frames) == 0:
                    paths.pop(key)

    log.message("Create separate paths:")

    pbar = log.pbar(len(paths), kind=pbar_name)
    # yield_single_paths requires a list of paths not a dictionary
    spaths = [sp for sp, nr in yield_single_paths(paths.values(), progress=True) if pbar.update(nr + 1) is None]

    pbar.finish()

    if options.sort_by_id:
        with log.fbm("Sort separate paths by resid"):
            spaths = sorted(spaths, key=lambda sp: sp.id)
    # apply smoothing?
    if options.apply_smoothing:
        log.message('Applying hard smoothing:')
        smooth = get_smooth_method(soptions)
        pbar = log.pbar(len(spaths), kind=pbar_name)
        for nr, sp in enumerate(spaths):
            sp.apply_smoothing(smooth)
            pbar.update(nr + 1)
        pbar.finish()
    if options.apply_soft_smoothing:
        log.message('Applying soft smoothing:')
        smooth = get_smooth_method(soptions)
        pbar = log.pbar(len(spaths), kind=pbar_name)
        for nr, sp in enumerate(spaths):
            sp.get_coords(smooth=smooth)
            pbar.update(nr + 1)
        pbar.finish()

    return {'paths': paths, 'spaths': spaths, 'options': options, 'soptions': soptions}


################################################################################

def stage_IV_run(config, options,
                 spaths=None,
                 **kwargs):
    coptions = config.get_cluster_options()
    # find coords of inlets, type and id
    inlet_coords = []
    inlet_type = []
    inlet_id = []
    inlet_spnr = []
    for nr, sp in enumerate(spaths):
        for inlet, itype in sp.coords_filo:
            inlet_coords.append(inlet.tolist())
            inlet_type.append(itype)
            inlet_id.append(sp.id)
            inlet_spnr.append(nr)
            # print nr,sp.id,inlet,itype

    # find nr of inlets
    nr_of_inlets = len(inlet_id)

    if nr_of_inlets > 0:
        with log.fbm("Performing clusterization"):
            assert coptions.method == 'dbscan', 'Unknown clusterization method %s.' % coptions.method
            # perform clusterization
            clusters = perform_clustering(np.array(inlet_coords), eps=float(coptions.eps),
                                          min_samples=int(coptions.min_samples))
            clusters += 1  # -1 is 0 now and it means unclustered
    else:
        log.message("No inlets found. Clusterization skipped.")
        clusters = np.array([])

    return {'inlet_coords': inlet_coords,
            'inlet_type': inlet_type,
            'inlet_id': inlet_id,
            'inlet_spnr': inlet_spnr,
            'clusters': clusters,
            'options': options,
            'coptions': coptions}


################################################################################

def stage_V_run(config, options,
                spaths=None,
                paths=None,
                inlet_id=None,
                clusters=None,
                inlet_type=None,
                **kwargs):
    # file handle?
    if options.save:
        fh = open(options.save, 'w')
        log.message('Using user provided file (%s).' % options.save)
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
    print >> fh, log.get_str_timestamp()

    ############
    if options.dump_config:
        print >> fh, asep()
        print >> fh, log.underline('Configuration file name: %s' % args.config_file)
        print >> fh, os.linesep.join(config.dump_config())

    ############
    print >> fh, asep()
    print >> fh, "Number of traceable residues:", len(paths)
    print >> fh, "Number of separate paths:", len(spaths)

    ############
    print >> fh, asep()
    print >> fh, "List of separate paths"
    header_template = " ".join(['%7s'] * 7)
    header = header_template % tuple("Nr ID Begin INP OBJ OUT End".split())
    print >> fh, log.thead(header)
    for nr, sp in enumerate(spaths):
        line = []
        for e in (nr, sp.id, sp.begins, len(sp.path_in), len(sp.path_object), len(sp.path_out), sp.ends):
            line += ["%7d" % e]
        print >> fh, " ".join(line)

    ############
    print >> fh, asep()
    print >> fh, "Separate paths lengths"
    header_template = " ".join(['%7s'] * 2 + ['%9s'] * 3)
    header = header_template % tuple("Nr ID INP OBJ OUT".split())
    print >> fh, log.thead(header)
    line_template = " ".join(['%7d'] * 2 + ['%9.1f'] * 3)
    for nr, sp in enumerate(spaths):
        line = [nr, sp.id]
        for e in sp.coords_in, sp.coords_object, sp.coords_out:
            if len(e) > 1:
                line += [sum(traces.diff(e))]
            else:
                line += [float('nan')]
        print >> fh, line_template % tuple(line)

    ############
    print >> fh, asep()
    print >> fh, "Separate paths average step lengths"
    header_template = " ".join(['%7s'] * 2 + ['%8s'] * 6)
    header = header_template % tuple("Nr ID INP INPstd OBJ OBJstd OUT OUTstd".split())
    print >> fh, log.thead(header)
    line_template = " ".join(['%7d'] * 2 + ['%8.2f', '%8.3f'] * 3)
    for nr, sp in enumerate(spaths):
        line = [nr, sp.id]
        for e in sp.coords_in, sp.coords_object, sp.coords_out:
            if len(e) > 1:
                line += [np.mean(traces.diff(e))]
                line += [np.std(traces.diff(e))]
            else:
                line += [float('nan')]
                line += [float('nan')]
        print >> fh, line_template % tuple(line)

    ############
    print >> fh, asep()
    print >> fh, "Number of inlets:", len(inlet_id)
    no_of_clusters = len(set([c for c in clusters if c != 0]))
    print >> fh, "Number of clusters:", no_of_clusters
    print >> fh, "Outliers:", {True: 'yes', False: 'no'}[0 in clusters]

    ############
    print >> fh, asep()
    print >> fh, "Clusters summary"
    clusters_list = list(set(clusters.tolist()))
    clusters_list.sort()
    header_template = " ".join(['%7s'] * 5)
    header = header_template % tuple("Nr Cluster Size INP OUT".split())
    print >> fh, log.thead(header)
    line_template = " ".join(['%7d'] * 5)
    for nr, c in enumerate(clusters_list):
        line = [nr, int(c), clusters.tolist().count(c)]
        # inlets types of current cluster
        its = [inlet_type[nr] for nr, cc in enumerate(clusters.tolist()) if cc == c]
        line.append(its.count(InletTypeCodes.inlet_in_code))
        line.append(its.count(InletTypeCodes.inlet_out_code))
        print >> fh, line_template % tuple(line)

    ############
    print >> fh, asep()
    print >> fh, "Separate paths inlets clusters"
    header_template = " ".join(['%7s'] * 2 + ['%7s'] * 2)
    header = header_template % tuple("Nr ID INP OUT".split())
    print >> fh, log.thead(header)
    line_template = " ".join(['%7d'] * 2 + ['%7s'] * 2)
    spaths_clust_type = []
    for nr, sp in enumerate(spaths):
        line = [nr, sp.id]
        ids = [nr for nr, iid in enumerate(inlet_id) if iid == sp.id]
        inp, out = None, None
        if len(ids) > 0:
            for iid in ids:
                if inlet_type[iid] == InletTypeCodes.inlet_in_code:
                    inp = int(clusters[iid])
                if inlet_type[iid] == InletTypeCodes.inlet_out_code:
                    out = int(clusters[iid])
        line.extend(map(str, [inp, out]))
        spaths_clust_type.append((inp, out))
        print >> fh, line_template % tuple(line)

    ############
    print >> fh, asep()
    print >> fh, "Clusters type summary"
    clusters_type_list = list(set(spaths_clust_type))
    clusters_type_list = sorted(sorted(clusters_type_list), key=lambda x: x[:1].count(None) + x.count(None))
    header_template = " ".join(['%7s'] * 3)
    header = header_template % tuple("Nr CType Size".split())
    print >> fh, log.thead(header)
    line_template = "%7d %7s %7d"

    def clusters_type_name(inp_cluster, out_cluster):
        name_of_cluster = ''
        if inp_cluster is None:
            name_of_cluster += 'N'
        else:
            name_of_cluster += str(inp_cluster)
        name_of_cluster += ':'
        if out_cluster is None:
            name_of_cluster += 'N'
        else:
            name_of_cluster += str(out_cluster)
        return name_of_cluster

    for nr, ct in enumerate(clusters_type_list):
        line = [nr]
        cts = clusters_type_name(*ct)
        line.append(cts)
        # inlets types of current cluster
        line.append(spaths_clust_type.count(ct))
        # calculate some statistics


        print >> fh, line_template % tuple(line)

    if options.save:
        fh.close()


################################################################################

def stage_VI_run(config, options,
                 spaths=None,
                 clusters=None,
                 inlet_coords=None,
                 **kwargs):
    from aqueduct.visual.pymol_connector import ConnectToPymol, SinglePathPlotter
    from aqueduct.visual.quickplot import ColorMapDistMap

    soptions = config.get_smooth_options()
    smooth = get_smooth_method(soptions)

    # start pymol
    with log.fbm("Starting PyMOL"):
        ConnectToPymol.init_pymol()
        spp = SinglePathPlotter()

    if options.all_paths_raw:
        with log.fbm("All raw paths"):
            if options.all_paths_split:
                spp.paths_trace(spaths, name='all_raw_in', plot_object=False, plot_out=False)
                spp.paths_trace(spaths, name='all_raw_obj', plot_in=False, plot_out=False)
                spp.paths_trace(spaths, name='all_raw_out', plot_in=False, plot_object=False)
            else:
                spp.paths_trace(spaths, name='all_raw')
    if options.all_paths_raw_io:
        spp.paths_inlets(spaths, name='all_raw_paths_io')

    if options.all_paths_smooth:
        with log.fbm("All smooth paths"):
            if options.all_paths_split:
                spp.paths_trace(spaths, name='all_smooth_in', plot_object=False, plot_out=False, smooth=smooth)
                spp.paths_trace(spaths, name='all_smooth_obj', plot_in=False, plot_out=False, smooth=smooth)
                spp.paths_trace(spaths, name='all_smooth_out', plot_in=False, plot_object=False, smooth=smooth)
            else:
                spp.paths_trace(spaths, name='all_smooth', smooth=smooth)
    if options.all_paths_smooth_io:
        spp.paths_inlets(spaths, name='all_smooth_paths_io', smooth=smooth)

    if options.paths_states:
        with log.fbm("Paths as states"):
            # as one object
            if options.paths_raw:
                [spp.paths_trace([sp], name='raw_paths', state=nr + 1) for nr, sp in enumerate(spaths)]
            if options.paths_smooth:
                [spp.paths_trace([sp], name='smooth_paths', smooth=smooth, state=nr + 1) for nr, sp in
                 enumerate(spaths)]
            if options.paths_raw_io:
                [spp.paths_inlets([sp], name='raw_paths_io', state=nr + 1) for nr, sp in enumerate(spaths)]
            if options.paths_smooth_io:
                [spp.paths_inlets([sp], name='smooth_paths_io', state=nr + 1, smooth=smooth) for nr, sp in
                 enumerate(spaths)]
    else:
        with log.fbm("Paths as separate objects"):
            if options.paths_raw:
                [spp.paths_trace([sp], name='raw_%d' % sp.id) for sp in spaths]
            if options.paths_smooth:
                [spp.paths_trace([sp], name='smooth_%d' % sp.id, smooth=smooth) for sp in spaths]
            if options.paths_raw_io:
                [spp.paths_inlets([sp], name='raw_paths_io_%d' % sp.id) for sp in spaths]
            if options.paths_smooth_io:
                [spp.paths_inlets([sp], name='smooth_paths_io_%d' % sp.id, smooth=smooth) for sp in spaths]

    if options.inlets_clusters:
        with log.fbm("Clusters"):
            # TODO: require stage V for that?
            no_of_clusters = len(set([c for c in clusters if c != -1]))
            cmap = ColorMapDistMap(name='Dark2', size=no_of_clusters)
            clusters_list = list(set(clusters.tolist()))
            clusters_list.sort()
            for c in clusters_list:
                # coords for current cluster
                ics = [inlet_coords[nr] for nr, cc in enumerate(clusters.tolist()) if cc == c]
                if c == -1:
                    c_name = 'none'
                else:
                    c_name = str(int(c))
                spp.scatter(ics, color=cmap(c), name="cluster_%s" % c_name)

    if options.show_molecule:
        with log.fbm("Molecule"):
            mda_ppr = mda.core.flags["permissive_pdb_reader"]
            mda.core.flags["permissive_pdb_reader"] = False
            pdb = TmpDumpWriterOfMDA()
            frames_to_show = range2int(options.show_molecule_frames)
            pdb.dump_frames(reader, frames=frames_to_show)
            ConnectToPymol.load_pdb('molecule', pdb.close())
            del pdb
            mda.core.flags["permissive_pdb_reader"] = mda_ppr


################################################################################


if __name__ == "__main__":
    ############################################################################
    # argument parsing
    import argparse

    parser = argparse.ArgumentParser(description="Valve, Aqueduct driver",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--dump-template-config", action="store_true", dest="dump_template_conf", required=False,
                        help="Dumps template config file. Supress all other output or actions.")

    parser.add_argument("-c", action="store", dest="config_file", required=False, help="Config file filename.")
    args = parser.parse_args()
    ############################################################################
    # special option for dumping template config
    config = ValveConfig()  # config template
    if args.dump_template_conf:
        import StringIO

        config_dump = StringIO.StringIO()
        config.save_config_stream(config_dump)
        print config_dump.getvalue()
        exit()

    ############################################################################
    # begin!

    valve_begin()
    valve_load_config(args.config_file, config)

    log.message("Optimal threads count: %d" % optimal_threads)
    # get global options
    goptions = config.get_global_options()
    pbar_name = goptions.pbar

    ############################################################################
    # STAGE 0

    reader = valve_read_trajectory(goptions.top, goptions.nc)

    ############################################################################
    # STAGE I
    max_frame = reader.number_of_frames - 1
    max_frame = 100
    result1 = valve_exec_stage(0, config, stage_I_run,
                               reader=reader,
                               max_frame=max_frame)

    ############################################################################
    # STAGE II
    result2 = valve_exec_stage(1, config, stage_II_run,
                               reader=reader,
                               max_frame=max_frame,
                               **result1)

    ############################################################################
    # STAGE III
    result3 = valve_exec_stage(2, config, stage_III_run,
                               **result2)

    ############################################################################
    # STAGE IV
    result4 = valve_exec_stage(3, config, stage_IV_run,
                               **result3)

    ############################################################################
    # STAGE V
    results = {}
    for result in (result2, result3, result4):
        results.update(result)

    result5 = valve_exec_stage(4, config, stage_V_run, no_io=True,
                               **results)

    ############################################################################
    # STAGE VI
    results = {}
    for result in (result3, result4):
        results.update(result)

    result6 = valve_exec_stage(5, config, stage_VI_run, no_io=True,
                               **results)

    ############################################################################
    # end!

    valve_end()
