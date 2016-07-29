#!/bin/env python2

"""
This is driver for aqueduct.
"""

import ConfigParser
import cPickle as pickle
import copy
import gzip
import multiprocessing as mp
import numpy as np
import os
import re
import shlex
import sys
from collections import namedtuple, OrderedDict
from functools import wraps
from itertools import izip_longest
from scipy.spatial.distance import cdist

import MDAnalysis as mda
import roman

from aqueduct import greetings as greetings_aqueduct
from aqueduct import version as aqueduct_version
from aqueduct import version_nice as aqueduct_version_nice
from aqueduct.geom import traces
from aqueduct.geom.cluster import PerformClustering, DBSCAN, AffinityPropagation, MeanShift, KMeans
from aqueduct.geom.convexhull import is_point_within_convexhull
from aqueduct.geom.master import create_master_spath
from aqueduct.geom.smooth import WindowSmooth, MaxStepSmooth, WindowOverMaxStepSmooth, ActiveWindowSmooth, \
    ActiveWindowOverMaxStepSmooth, DistanceWindowSmooth, DistanceWindowOverMaxStepSmooth
from aqueduct.traj.dumps import TmpDumpWriterOfMDA
from aqueduct.traj.inlets import Inlets, InletTypeCodes
from aqueduct.traj.paths import GenericPaths, yield_single_paths
from aqueduct.traj.reader import ReadViaMDA
from aqueduct.traj.selections import CompactSelectionMDA, SelectionMDA
from aqueduct.utils import log
from aqueduct.utils.helpers import range2int, Auto, what2what, lind

# TODO: Move it to separate module
cpu_count = mp.cpu_count()

# global optimal_threads
optimal_threads = None


def version():
    return 0, 8, 1


def version_nice():
    return '.'.join(map(str, version()))


__mail__ = 'info@aquaduct.pl'
__version__ = version_nice()


###############################################################################
# configuration file helpers

class ConfigSpecialNames:
    special_names_dict = {'none': None,
                          'null': None,
                          'true': True,
                          'false': False,
                          'yes': True,
                          'no': False,
                          'auto': Auto}

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
        return 'execute dump'.split()

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
    def recluster_name():
        return 'reclusterization'

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

    def get_recluster_options(self):
        section = self.recluster_name()
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
        #snr = 0 # stage number

        config = ConfigParser.RawConfigParser()

        def common(section):
            for setting in self.common_config_names():
                value = None
                if setting == 'execute':
                    value = 'runonce'
                elif setting == 'dump':
                    value = '%d_%s_data.dump' % (snr+1,section)
                if value is None:
                    config.set(section, setting)
                else:
                    config.set(section, setting, value=value)

        def common_traj_data(section):
            for setting in self.common_traj_data_config_names():
                config.set(section, setting)

        ################
        # global settings
        section = self.global_name()
        config.add_section(section)

        # top - top file name
        # nc - netcdf file name
        config.set(section, 'top') # topology
        config.set(section, 'trj') # trajectory

        #config.set(section, 'pbar', 'simple')


        ################
        snr = 0 # stage number
        # stage I
        # find traceable residues
        section = self.stage_names(snr)
        config.add_section(section)

        common(section)
        common_traj_data(section)

        ################
        snr+=1
        # stage II
        # find raw paths
        section = self.stage_names(snr)
        config.add_section(section)

        common(section)
        common_traj_data(section)

        config.set(section, 'clear_in_object_info', 'False')

        ################
        snr+=1
        # stage III
        # create separate frames
        section = self.stage_names(snr)
        config.add_section(section)

        common(section)

        config.set(section, 'discard_empty_paths', 'True')
        config.set(section, 'sort_by_id', 'True')
        config.set(section, 'apply_smoothing', 'False')
        config.set(section, 'apply_soft_smoothing', 'True')
        config.set(section, 'discard_short_paths', '1')

        ################
        snr+=1
        # stage IV
        # inlets clusterisation
        section = self.stage_names(snr)
        config.add_section(section)

        common(section)

        config.set(section, 'recluster_outliers', 'False')
        config.set(section, 'detect_outliers', 'False')
        config.set(section, 'singletons_outliers', 'False')

        ################
        # smooth
        section = self.smooth_name()
        config.add_section(section)
        config.set(section, 'method', 'window')

        ################
        # clusterization
        section = self.cluster_name()
        config.add_section(section)
        config.set(section, 'method', 'meanshift')
        config.set(section, 'cluster_all', 'False')
        config.set(section, 'bandwidth', 'Auto')

        ################
        # reclusterization
        section = self.recluster_name()
        config.add_section(section)
        config.set(section, 'method', 'dbscan')
        config.set(section, 'eps', '5.0')
        config.set(section, 'min_samples', '3')

        ################
        snr+=1
        # stage V
        # analysis
        section = self.stage_names(snr)
        config.add_section(section)
        common(section)
        config.remove_option(section, 'dump')
        config.set(section, 'save', '%d_%s_results.txt' % (snr + 1, section))

        config.set(section, 'dump_config', 'True')

        ################
        snr+=1
        # stage VI
        # visualize
        section = self.stage_names(snr)
        config.add_section(section)
        common(section)
        config.remove_option(section, 'dump')
        config.set(section, 'save', '%d_%s_results.pse' % (snr+1,section))

        config.set(section, 'simply_smooths', 0.05236)

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

        config.set(section, 'ctypes_raw', 'True')
        config.set(section, 'ctypes_smooth', 'True')

        # visualize clusters
        config.set(section, 'inlets_clusters', 'True')

        # show protein
        config.set(section, 'show_molecule', 'None')
        config.set(section, 'show_molecule_frames', '0')
        config.set(section, 'show_chull', 'None')
        config.set(section, 'show_chull_frames', '0')

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
                  [self.get_stage_options(stage) for stage in range(6)] + \
                  [self.get_cluster_options(), self.get_recluster_options(),
                   self.get_smooth_options()]
        names = [self.global_name()] + \
                [self.stage_names(stage) for stage in range(6)] + \
                [self.cluster_name(), self.recluster_name(),
                 self.smooth_name()]

        for o, n in zip(options, names):
            output.append('[%s]' % n)
            for k in o._asdict().keys():
                v = o._asdict()[k]
                if v is Auto:  # FIXME: do something with Auto class!
                    v = str(Auto())
                else:
                    v = str(v)
                output.append('%s = %s' % (k, v))

        return output


################################################################################
# reader helper class

class TrajectoryReader(object):
    def __init__(self, top, trj):
        assert isinstance(top, (str, unicode)), "Topology file name missing, %s given instead" % str(top)
        assert isinstance(trj, (str, unicode)), "Trajectory file(s) name(s) missing, %s given instead" % str(trj)
        self.top = top
        self.trj = shlex.split(trj)

    def get(self):
        # assume it is a Amber NetCDF
        # TODO: check if it is DCD and do something?
        # TODO: move it to another class, ReaderHelper for instance.

        return ReadViaMDA(self.top, self.trj)
        # return ReadAmberNetCDFviaMDA(self.top, self.trj)

    @property
    def max_frame(self):
        with self.get() as tmp_reader:
            return tmp_reader.number_of_frames - 1  # returns 0-based value


def rebuild_selection(selection, reader):
    return CompactSelectionMDA(selection).toSelectionMDA(reader)


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
        is_res_in_scope = CHullCheck_exec(chull, res_coords, threads=optimal_threads)
    else:
        if res.unique_resids_number() == 0:
            return []
        res_in_scope_uids = scope.unique_resids(ikwid=True)
        # res_in_scope_uids = traj_reader.parse_selection(options.scope).unique_resids(ikwid=True)
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
        log.error('Loaded data has %s major version lower then the application.' % what)
    if current[0] < loaded[0]:
        log.error('Loaded data has %s major version higher then the application.' % what)
    if current[0] != loaded[0]:
        log.error('Possible problems with API compilance.')
    if current[1] > loaded[1]:
        log.warning('Loaded data has %s minor version lower then the application.' % what)
    if current[1] < loaded[1]:
        log.warning('Loaded data has %s minor version higher then the application.' % what)
    if current[1] != loaded[1]:
        log.warning('Possible problems with API compilance.')


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
                with reader.get() as traj_reader:
                    value = value.toSelectionMDA(traj_reader)
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


def get_clustering_method(coptions):
    assert coptions.method in ['dbscan', 'affprop', 'meanshift',
                               'kmeans'], 'Unknown clusterization method %s.' % coptions.method

    opts = {}

    def dbscan_opts():
        if 'eps' in coptions._asdict():
            opts.update({'eps': float(coptions.eps)})
        if 'min_samples' in coptions._asdict():
            opts.update({'min_samples': int(coptions.min_samples)})

    def affprop_opts():
        pass

    def kmeans_opts():
        if 'n_clusters' in coptions._asdict():
            opts.update({'n_clusters': int(coptions.n_clusters)})

    def meanshift_opts():
        if 'cluster_all' in coptions._asdict():
            opts.update({'cluster_all': bool(coptions.cluster_all)})
        if 'bandwidth' in coptions._asdict():
            if coptions.bandwidth in (Auto, None):
                opts.update({'bandwidth': coptions.bandwidth})
            else:
                opts.update({'bandwidth': float(coptions.bandwidth)})

    if coptions.method == 'dbscan':
        dbscan_opts()
        method = DBSCAN
    elif coptions.method == 'affprop':
        affprop_opts()
        method = AffinityPropagation
    elif coptions.method == 'kmeans':
        kmeans_opts()
        method = KMeans
    elif coptions.method == 'meanshift':
        meanshift_opts()
        method = MeanShift

    return PerformClustering(method, **opts)
    # return lambda X: perform_clustering(X, method, **opts)


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
    assert os.path.isfile(filename), "Config file %s does not exist." % filename
    with log.fbm('Load configuration file'):
        config.load_config(filename)


def valve_read_trajectory(top, traj):
    with log.fbm('Read trajectory'):
        return TrajectoryReader(top, traj)
        # read trajectory
        # traj_list = shlex.split(traj)
        # return ReadAmberNetCDFviaMDA(top, traj_list)
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

    # TODO: consder to create traj_reader object here instead of doing it in stage_run or in load...
    # execute?
    can_be_loaded = False
    if (not no_io) and options.dump:
        if os.path.isfile(options.dump) or os.path.islink(options.dump):
            can_be_loaded = True
    if options.execute in ['run'] or (options.execute in ['runonce'] and not can_be_loaded):
        result = stage_run(config, options, reader=reader, **kwargs)
        if not no_io:
            ###########
            # S A V E #
            ###########
            save_stage_dump(options.dump, **result)
    elif options.execute in ['skip'] or (options.execute in ['runonce'] and can_be_loaded):
        if not no_io:
            ###########
            # L O A D #
            ###########
            if options.dump:
                result = load_stage_dump(options.dump, reader=reader)
    else:
        raise NotImplementedError('exec mode %s not implemented' % options.execute)
    # remove options stuff
    if not no_io:
        if result is not None:
            return dict(((key, val) for key, val in result.iteritems() if 'options' not in key))


################################################################################
# stages run

# traceable_residues
def stage_I_run(config, options,
                reader=None,
                max_frame=None,
                **kwargs):
    log.message("Loop over frames - search of residues in object:")
    pbar = log.pbar(max_frame)

    # create pool of workers - mapping function
    map_fun = map
    if optimal_threads > 1:
        pool = mp.Pool(optimal_threads)
        map_fun = pool.map

    with reader.get() as traj_reader:

        scope = traj_reader.parse_selection(options.scope)

        # create some containers
        res_ids_in_object_over_frames = {}
        all_res = None

        for frame in traj_reader.iterate_over_frames():
            if frame > max_frame:
                break
            # current res selection
            res = traj_reader.parse_selection(options.object)

            # check if res are in scope
            if options.scope_convexhull:
                res_coords = list(res.center_of_mass_of_residues())
                chull = scope.get_convexhull_of_atom_positions()
                is_res_in_scope = map_fun(is_point_within_convexhull, izip_longest(res_coords, [], fillvalue=chull))
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
    # destroy pool of workers
    if optimal_threads > 1:
        pool.close()
        pool.join()
        del pool

    pbar.finish()

    if all_res is None:
        raise ValueError("No traceable residues was found.")

    log.message("Number of residues to trace: %d" % all_res.unique_resids_number())

    return {'all_res': all_res,
            'res_ids_in_object_over_frames': res_ids_in_object_over_frames,
            'options': options}


################################################################################

# raw_paths
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

    with reader.get() as traj_reader:

        scope = traj_reader.parse_selection(options.scope)

        with log.fbm("Rebuild treceable residues with current trajectory"):
            all_res = rebuild_selection(all_res, traj_reader)

        log.message("Trajectory scan:")
        pbar = log.pbar(max_frame)

        # create pool of workers - mapping function
        map_fun = map
        if optimal_threads > 1:
            pool = mp.Pool(optimal_threads)
            map_fun = pool.map

        for frame in traj_reader.iterate_over_frames():
            if frame > max_frame:
                break

            all_res_coords = list(all_res.center_of_mass_of_residues())  # this uses iterate over residues

            # check if is res are in scope
            # is_res_in_scope = check_res_in_scope(options, scope, all_res, all_res_coords)
            # check is res are in scope
            if options.scope_convexhull:
                chull = scope.get_convexhull_of_atom_positions()
                is_res_in_scope = map_fun(is_point_within_convexhull, izip_longest(all_res_coords, [], fillvalue=chull))
            else:
                is_res_in_scope = check_res_in_scope(options, scope, res, None)

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
                        res = traj_reader.parse_selection(options.object)
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

    # destroy pool of workers
    if optimal_threads > 1:
        pool.close()
        pool.join()
        del pool

    pbar.finish()

    return {'all_res': all_res, 'paths': paths, 'options': options}


################################################################################

# separate_paths
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
    pbar = log.pbar(len(paths))
    # yield_single_paths requires a list of paths not a dictionary
    spaths = [sp for sp, nr in yield_single_paths(paths.values(), progress=True) if pbar.update(nr + 1) is None]
    pbar.finish()

    if options.discard_short_paths > 0:
        shorter_then = int(options.discard_short_paths)
        with log.fbm("Discard paths shorter then %d" % shorter_then):
            spaths = [sp for sp in spaths if sp.size > shorter_then]

    if options.sort_by_id:
        with log.fbm("Sort separate paths by resid"):
            spaths = sorted(spaths, key=lambda sp: (sp.id.id, sp.id.nr))
    # apply smoothing?
    if options.apply_smoothing or options.apply_soft_smoothing:
        smooth = get_smooth_method(soptions)
    if options.apply_smoothing:
        log.message('Applying hard smoothing:')
        pbar = log.pbar(len(spaths))
        for nr, sp in enumerate(spaths):
            sp.apply_smoothing(smooth)
            pbar.update(nr + 1)
        pbar.finish()
    if options.apply_soft_smoothing:
        log.message('Applying soft smoothing:')
        pbar = log.pbar(len(spaths))
        for nr, sp in enumerate(spaths):
            sp.get_coords(smooth=smooth)
            pbar.update(nr + 1)
        pbar.finish()

    return {'paths': paths, 'spaths': spaths, 'options': options, 'soptions': soptions}


################################################################################

# inlets_clusterization
def stage_IV_run(config, options,
                 spaths=None,
                 **kwargs):
    coptions = config.get_cluster_options()
    rcoptions = config.get_recluster_options()
    soptions = config.get_smooth_options()

    # new style clustering
    with log.fbm("Create inlets"):
        inls = Inlets(spaths)

    def noo():
        # returns number of outliers
        if 0 in inls.clusters_list:
            return inls.clusters.count(0)
        return 0

    if inls.size > 0:
        # ***** CLUSTERIZATION *****
        with log.fbm("Performing clusterization"):
            clustering_function = get_clustering_method(coptions)
            inls.perform_clustering(clustering_function)
        log.message('Number of clusters detected so far: %d' % len(inls.clusters_list))
        # ***** OUTLIERS DETECTION *****
        if options.detect_outliers:
            log.message('Number of outliers so far: %d' % noo())
            with log.fbm("Detecting outliers"):
                if options.detect_outliers is not Auto:
                    threshold = float(options.detect_outliers)
                clusters = inls.clusters
                for cluster, center, std in zip(inls.clusters_list,
                                                inls.clusters_centers,
                                                inls.clusters_std):
                    if cluster == 0:
                        continue
                    inls_lim = inls.lim2clusters(cluster)
                    dmat = cdist(np.matrix(center), inls_lim.coords, metric='euclidean').flatten()
                    # Auto procedure
                    if options.detect_outliers is Auto:
                        threshold = std * 4  # FIXME: magic constant!
                    # defined threshold procedure
                    for nr, (d, ids) in enumerate(zip(dmat, inls_lim.inlets_ids)):
                        # print d, threshold
                        if d > threshold:
                            clusters[ids] = 0
                inls.add_cluster_annotations(clusters)
        log.message('Number of clusters detected so far: %d' % len(inls.clusters_list))
        log.message('Number of outliers: %d' % noo())
        # ***** RECLUSTERIZATION *****
        if options.recluster_outliers:
            with log.fbm("Performing reclusterization of outliers"):
                clustering_function = get_clustering_method(rcoptions)
                # perform reclusterization
                inls.recluster_outliers(clustering_function)
            log.message('Number of clusters detected so far: %d' % len(inls.clusters_list))
            log.message('Number of outliers: %d' % noo())
        # ***** SINGLETONS REMOVAL *****
        if options.singletons_outliers:
            with log.fbm("Removing clusters of size %d" % int(options.singletons_outliers)):
                inls.small_clusters_to_outliers(int(options.singletons_outliers))
            log.message('Number of clusters detected so far: %d' % len(inls.clusters_list))
            log.message('Number of outliers: %d' % noo())

        with log.fbm("Calculating cluster types"):
            ctypes = inls.spaths2ctypes(spaths)

        # now, there is something to do with ctypes!
        # we can create master paths!

        log.message("Creating master paths for cluster types:")

        smooth = get_smooth_method(soptions)
        master_paths = {}
        master_paths_smooth = {}
        ctypes_generic = [ct.generic for ct in ctypes]
        ctypes_generic_list = sorted(list(set(ctypes_generic)))

        pbar = log.pbar(len(ctypes_generic_list) * 2)

        for nr, ct in enumerate(ctypes_generic_list):
            sps = lind(spaths, what2what(ctypes_generic, [ct]))
            # print len(sps),ct
            master_paths.update({ct: create_master_spath(sps, resid=nr, ctype=ct, heartbeat=pbar.heartbeat)})
            pbar.update(nr * 2)
            master_paths_smooth.update(
                {ct: create_master_spath(sps, resid=nr, ctype=ct, smooth=smooth, heartbeat=pbar.heartbeat)})
            pbar.update(nr * 2 + 1)
        pbar.finish()

    else:
        log.message("No inlets found. Clusterization skipped.")
        # make empty results
        ctypes = inls.spaths2ctypes(spaths)
        master_paths = {}
        master_paths_smooth = {}

    return {'inls': inls,
            'ctypes': ctypes,
            'master_paths': master_paths,
            'master_paths_smooth': master_paths_smooth}


################################################################################

def make_line(template, line):
    return (' '.join(template)) % tuple(line)


def make_header_template(line_template):
    header_template = []
    col_re = re.compile('[0-9]+')
    for l in line_template:
        header_template.append('%{0}s'.format(col_re.findall(l)[0]))
    return header_template


def nr_header():
    return ['Nr'], ['%7d']


def get_header_line_and_line_template(header_line_and_line_template, head_nr=False):
    header, line_template = header_line_and_line_template

    header_template = make_header_template(line_template)
    # head_nr? only to header
    if head_nr:
        nrh, nrlt = nr_header()
        header_template = make_header_template(nrlt + line_template)
        header = nrh + header
    header_line = make_line(header_template, header)

    return header_line, line_template


def spath_id_header():
    return ['ID'], ['%9s']


def add_path_id_head(gen):
    sph, splt = spath_id_header()

    @wraps(gen)
    def patched(*args, **kwargs):
        add_id = True
        if 'add_id' in kwargs:
            add_id = kwargs.pop('add_id')
        h, lt = gen(*args, **kwargs)
        if add_id:
            return sph + h, splt + lt
        return h, lt

    return patched


def add_path_id(gen):
    @wraps(gen)
    def patched(spath, *args, **kwargs):
        add_id = True
        if 'add_id' in kwargs:
            add_id = kwargs.pop('add_id')
        line = gen(spath, *args, **kwargs)
        if add_id:
            line = [spath.id] + line
        return line

    return patched


def size_header():
    return ['Size'], ['%7d']


def add_size_head(gen):
    sph, splt = size_header()

    @wraps(gen)
    def patched(*args, **kwargs):
        h, lt = gen(*args, **kwargs)
        return sph + h, splt + lt

    return patched


def add_size(gen):
    @wraps(gen)
    def patched(spaths, add_size=True, *args, **kwargs):
        line = gen(spaths, *args, **kwargs)
        if add_size:
            line = [len(spaths)] + line
        return line

    return patched


def cluster_id_header():
    return ['Cluster'], ['%7d']


def add_cluster_id_head(gen):
    sph, splt = cluster_id_header()

    @wraps(gen)
    def patched(*args, **kwargs):
        h, lt = gen(*args, **kwargs)
        return sph + h, splt + lt

    return patched


def add_cluster_id(gen):
    @wraps(gen)
    def patched(cluster, something, add_id=True, *args, **kwargs):
        line = gen(cluster, something, *args, **kwargs)
        if add_id:
            line = [int(cluster)] + line
        return line

    return patched


def ctype_id_header():
    return ['CType'], ['%7s']


def add_ctype_id_head(gen):
    sph, splt = ctype_id_header()

    @wraps(gen)
    def patched(*args, **kwargs):
        h, lt = gen(*args, **kwargs)
        return sph + h, splt + lt

    return patched


def add_ctype_id(gen):
    @wraps(gen)
    def patched(ctype, something, add_id=True, *args, **kwargs):
        line = gen(ctype, something, *args, **kwargs)
        if add_id:
            line = [str(ctype)] + line
        return line

    return patched


class PrintAnalysis(object):
    nr_template = '%7d '

    def __init__(self, fileoption):
        self.output2stderr = False
        if fileoption:
            self.filehandle = open(fileoption, 'w')
            # self.output2stderr = True
        else:
            self.filehandle = sys.stdout

    def __call__(self, info2print, nr=None):
        if nr is not None:
            info2print = (self.nr_template % nr) + info2print
        if self.output2stderr:
            log.message(info2print)
        print >> self.filehandle, info2print

    def sep(self):
        self(asep())

    def thead(self, info2print):
        self(log.thead(info2print))

    def under(self, info2print):
        self(log.underline(info2print))


################################################################################

@add_path_id_head
def spath_basic_info_header():
    header = 'BeginF InpF ObjF OutF EndF'.split()
    line_template = ['%7d'] * len(header)
    return header, line_template


@add_path_id
def spath_basic_info(spath):
    line = []
    line.append(spath.begins)
    line.extend(map(len, (spath.path_in, spath.path_object, spath.path_out)))
    line.append(spath.ends)
    return line


################

@add_path_id_head
def spath_lenght_total_info_header():
    header = 'InpL ObjL OutL'.split()
    line_template = ['%9.1f'] * len(header)
    return header, line_template


@add_path_id
def spath_lenght_total_info(spath):
    line = []
    for t in traces.midpoints(spath.coords):
        if len(t) > 0:
            line.append(traces.length_step_std(t)[0])
        else:
            line.append(float('nan'))
    return line


################

@add_path_id_head
def spath_steps_info_header():
    header = 'InpS InpStdS ObjS ObjStdS OutS OutStdS'.split()
    line_template = ['%8.2f', '%8.3f'] * (len(header) / 2)
    return header, line_template


@add_path_id
def spath_steps_info(spath):
    line = []
    for t in traces.midpoints(spath.coords):
        if len(t) > 0:
            line.extend(traces.length_step_std(t)[1:])
        else:
            line.extend([float('nan'), float('nan')])
    return line


################

@add_path_id_head
def spath_ctype_header():
    header, line_template = ctype_id_header()
    return header, line_template


@add_path_id
def spath_ctype(spath, ctype=None):
    line = [str(ctype)]
    return line


################

@add_path_id_head
def spath_full_info_header():
    header = []
    line_template = []
    for h, lt in (spath_basic_info_header(add_id=False),
                  spath_lenght_total_info_header(add_id=False),
                  spath_steps_info_header(add_id=False),
                  spath_ctype_header(add_id=False)):
        header += h
        line_template += lt
    return header, line_template


@add_path_id
def spath_full_info(spath, ctype=None):
    line = []
    for l in (spath_basic_info(spath, add_id=False),
              spath_lenght_total_info(spath, add_id=False),
              spath_steps_info(spath, add_id=False),
              spath_ctype(spath, ctype=ctype, add_id=False)):
        line += l
    return line


################################################################################

@add_size_head
def spaths_lenght_total_header():
    header = 'Inp InpStd Obj ObjStd Out OutStd'.split()
    line_template = ['%8.1f', '%8.2f'] * (len(header) / 2)
    return header, line_template


@add_size
def spaths_length_total(spaths):
    line = []
    d4s = []
    for sp in spaths:
        d4s.append(spath_lenght_total_info(sp, add_id=False))
    d4s = np.array(d4s)
    line.extend(np.mean(d4s, 0))
    line.extend(np.std(d4s, 0))
    return [line[0], line[3], line[1], line[4], line[2], line[5]]


################################################################################


@add_cluster_id_head
def clusters_inlets_header():
    header = 'Size SInp IInp IOut SOut'.split()
    header = 'Size INCOMING OUTGOING'.split()
    line_template = ['%7d'] + ['%8d'] * (len(header) - 1)
    return header, line_template


@add_cluster_id
def clusters_inlets(cluster, inlets):
    line = [inlets.size]
    line.append(inlets.lim2types([InletTypeCodes.surface_incoming]).size)
    # line.append(inlets.lim2types([InletTypeCodes.internal_incoming]).size)
    # line.append(inlets.lim2types([InletTypeCodes.internal_outgoing]).size)
    line.append(inlets.lim2types([InletTypeCodes.surface_outgoing]).size)
    return line


@add_ctype_id_head
def ctypes_spaths_info_header():
    header, line_template = spaths_lenght_total_header()
    return header, line_template


@add_ctype_id
def ctypes_spaths_info(ctype, spaths):
    line = []
    line += spaths_length_total(spaths)
    return line


################################################################################


# analysis
def stage_V_run(config, options,
                spaths=None,
                paths=None,
                inls=None,
                ctypes=None,
                **kwargs):
    # file handle?
    pa = PrintAnalysis(options.save)
    if options.save:
        log.message('Using user provided file (%s).' % options.save)
        log.message(sep())
        log.message('')
    else:
        log.message('Using standard output.')
        log.message(sep())
        log.message('')

    ############
    pa.sep()
    pa('Aqueduct analysis')
    pa(log.get_str_timestamp())

    ############
    if options.dump_config:
        pa.sep()
        pa.under('Configuration file name: %s' % args.config_file)
        pa(os.linesep.join(config.dump_config()))

    ############
    pa.sep()
    pa("Number of traceable residues: %d" % len(paths))
    pa("Number of separate paths: %d" % len(spaths))

    ############
    pa.sep()
    pa("Number of inlets: %d" % inls.size)
    no_of_clusters = len(inls.clusters_list) - {True: 1, False: 0}[0 in inls.clusters_list]  # minus outliers, if any
    pa("Number of clusters: %d" % no_of_clusters)
    pa("Outliers: %s" % ({True: 'yes', False: 'no'}[0 in inls.clusters_list]))

    ############
    pa.sep()
    pa("Clusters summary - inlets")
    header_line, line_template = get_header_line_and_line_template(clusters_inlets_header(), head_nr=True)
    pa.thead(header_line)

    for nr, cl in enumerate(inls.clusters_list):
        inls_lim = inls.lim2clusters(cl)
        pa(make_line(line_template, clusters_inlets(cl, inls_lim)), nr=nr)

    ############
    pa.sep()
    pa("Separate paths clusters types summary - mean lengths of paths")
    header_line, line_template = get_header_line_and_line_template(ctypes_spaths_info_header(), head_nr=True)
    pa.thead(header_line)

    ctypes_generic = [ct.generic for ct in ctypes]
    ctypes_generic_list = sorted(list(set(ctypes_generic)))

    # sorted by ctype
    ctypes_size = []
    for nr, ct in enumerate(ctypes_generic_list):
        sps = lind(spaths, what2what(ctypes_generic, [ct]))
        ctypes_size.append(len(sps))
        # pa(make_line(line_template, ctypes_spaths_info(ct, sps)), nr=nr)

    # sorted by sizes:
    ctypes_generic_list = sorted(ctypes_generic_list, key=lambda ctyp: ctypes_size[ctypes_generic_list.index(ctyp)],
                                 reverse=True)
    for nr, ct in enumerate(ctypes_generic_list):
        sps = lind(spaths, what2what(ctypes_generic, [ct]))
        ctypes_size.append(len(sps))
        pa(make_line(line_template, ctypes_spaths_info(ct, sps)), nr=nr)

    ############
    pa.sep()
    pa("List of separate paths and properties")
    header_line, line_template = get_header_line_and_line_template(spath_full_info_header(), head_nr=True)
    pa.thead(header_line)
    for nr, (sp, ctype) in enumerate(izip_longest(spaths, ctypes, fillvalue=None)):
        if ctype is not None:
            ctype = ctype.generic
        pa(make_line(line_template, spath_full_info(sp, ctype=ctype)), nr=nr)


################################################################################

# visualize


def plot_spaths_traces(spaths, spp=None, name=None, split=False, states=False, separate=False, smooth=None):
    if states or separate:
        spaths_iter = spaths
    else:
        spaths_iter = [spaths]
    state = None
    name_separate = ''
    for nr, sp in enumerate(spaths_iter):
        if states:
            state = nr + 1
        if separate:
            name_separate = '_%d' % (nr + 1)
        if states or separate:
            sp = [sp]
        if split:
            spp.paths_trace(sp, name=name + '_in' + name_separate, plot_object=False, plot_out=False, state=state,
                            smooth=smooth)
            spp.paths_trace(sp, name=name + '_obj' + name_separate, plot_in=False, plot_out=False, state=state,
                            smooth=smooth)
            spp.paths_trace(sp, name=name + '_out' + name_separate, plot_in=False, plot_object=False, state=state,
                            smooth=smooth)
        else:
            spp.paths_trace(sp, name=name + name_separate, state=state, smooth=smooth)


def plot_spaths_inlets(spaths, spp=None, name=None, states=False, separate=False, smooth=None):
    if states or separate:
        spaths_iter = spaths
    else:
        spaths_iter = [spaths]
    state = None
    name_separate = ''
    for nr, sp in enumerate(spaths_iter):
        if states:
            state = nr + 1
        if separate:
            name_separate = '_%d' % (nr + 1)
        if states or separate:
            sp = [sp]
        spp.paths_inlets(sp, name=name + name_separate, state=state, smooth=smooth)


def stage_VI_run(config, options,
                 reader=None,
                 spaths=None,
                 inls=None,
                 ctypes=None,
                 master_paths=None,
                 master_paths_smooth=None,
                 **kwargs):
    from aqueduct.visual.pymol_connector import ConnectToPymol, SinglePathPlotter
    from aqueduct.visual.pymol_connector import cmd as pymol_cmd
    from aqueduct.visual.quickplot import ColorMapDistMap

    soptions = config.get_smooth_options()
    smooth = get_smooth_method(soptions)

    # start pymol
    with log.fbm("Starting PyMOL"):
        # TODO: ConnectToPymol is used to initialize PyMol and to load molecule
        # TODO: SinglePathPlotter is used to put paths to PyMol
        # TODO: Both can be bassically changed in such a way that appropriate pdb files
        # TODO: would be generated and a companion script that would load them to PyMol
        pymol_connector = ConnectToPymol()
        #pymol_connector.init_pymol()
        pymol_connector.init_script('test_.py')

        if options.simply_smooths:
            spp = SinglePathPlotter(pymol_connector,linearize=float(options.simply_smooths))
        else:
            spp = SinglePathPlotter(pymol_connector,linearize=None)

    ctypes_generic = [ct.generic for ct in ctypes]
    ctypes_generic_list = sorted(list(set(ctypes_generic)))

    if options.show_molecule:
        with log.fbm("Molecule"):
            with reader.get() as traj_reader:
                mda_ppr = mda.core.flags["permissive_pdb_reader"]
                mda.core.flags["permissive_pdb_reader"] = False
                pdb = TmpDumpWriterOfMDA()
                frames_to_show = range2int(options.show_molecule_frames)
                pdb.dump_frames(traj_reader, frames=frames_to_show, selection=options.show_molecule)
                pymol_connector.load_pdb('molecule', pdb.close())
                del pdb
                mda.core.flags["permissive_pdb_reader"] = mda_ppr
                # it would be nice to plot convexhull
    if options.show_chull:
        with log.fbm("Convexhull"):
            with reader.get() as traj_reader:
                options_stageII = config.get_stage_options(1)
                if options_stageII.scope_convexhull:
                    scope = traj_reader.parse_selection(options_stageII.scope)
                    frames_to_show = range2int(options.show_chull_frames)
                    for frame in frames_to_show:
                        traj_reader.set_current_frame(frame)
                        chull = scope.get_convexhull_of_atom_positions()
                        spp.convexhull(chull,state=frame+1)


    if options.inlets_clusters:
        with log.fbm("Clusters"):
            # TODO: require stage V for that?
            no_of_clusters = len(inls.clusters_list)  # total, including outliers
            cmap = ColorMapDistMap(size=no_of_clusters)
            for c in inls.clusters_list:
                # coords for current cluster
                ics = inls.lim2clusters(c).coords
                if c == 0:
                    c_name = 'out'
                else:
                    c_name = str(int(c))
                spp.scatter(ics, color=cmap(c), name="cluster_%s" % c_name)

    if options.ctypes_raw:
        with log.fbm("CTypes raw"):
            for nr, ct in enumerate(ctypes_generic_list):
                log.message(str(ct),cont=True)
                sps = lind(spaths, what2what(ctypes_generic, [ct]))
                plot_spaths_traces(sps, name=str(ct) + '_raw', split=False, spp=spp)
                if ct in master_paths:
                    plot_spaths_traces([master_paths[ct]], name=str(ct) + '_raw_master', split=False, spp=spp)

    if options.ctypes_smooth:
        with log.fbm("CTypes smooth"):
            for nr, ct in enumerate(ctypes_generic_list):
                log.message(str(ct),cont=True)
                sps = lind(spaths, what2what(ctypes_generic, [ct]))
                plot_spaths_traces(sps, name=str(ct) + '_smooth', split=False, spp=spp, smooth=smooth)
                if ct in master_paths_smooth:
                    plot_spaths_traces([master_paths_smooth[ct]], name=str(ct) + '_smooth_master', split=False, spp=spp,
                                       smooth=lambda anything: anything)
                if ct in master_paths:
                    plot_spaths_traces([master_paths[ct]], name=str(ct) + '_raw_master_smooth', split=False, spp=spp,
                                       smooth=smooth)

    if options.all_paths_raw:
        with log.fbm("All raw paths"):
            plot_spaths_traces(spaths, name='all_raw', split=options.all_paths_split, spp=spp)
    if options.all_paths_raw_io:
        with log.fbm("All raw paths io"):
            plot_spaths_inlets(spaths, name='all_raw_paths_io', spp=spp)

    if options.all_paths_smooth:
        with log.fbm("All smooth paths"):
            plot_spaths_traces(spaths, name='all_smooth', split=options.all_paths_split, spp=spp, smooth=smooth)
    if options.all_paths_smooth_io:
        with log.fbm("All smooth paths io"):
            plot_spaths_inlets(spaths, name='all_smooth_paths_io', spp=spp)

    with log.fbm("Paths as states"):
        if options.paths_raw:
            log.message("raw", cont=True)
            plot_spaths_traces(spaths, name='raw_paths', states=options.paths_states, separate=not options.paths_states,
                               spp=spp)
        if options.paths_smooth:
            log.message("smooth", cont=True)
            plot_spaths_traces(spaths, name='smooth_paths', states=options.paths_states,
                               separate=not options.paths_states, smooth=smooth, spp=spp)
        if options.paths_raw_io:
            log.message("raw_io", cont=True)
            plot_spaths_inlets(spaths, name='raw_paths_io', states=options.paths_states,
                               separate=not options.paths_states, spp=spp)
        if options.paths_smooth_io:
            log.message("smooth_io", cont=True)
            plot_spaths_inlets(spaths, name='smooth_paths_io', states=options.paths_states,
                               separate=not options.paths_states, smooth=smooth, spp=spp)

    if options.show_molecule:
        pymol_connector.orient_on('molecule')

    if options.save:
        with log.fbm("Saving session (%s)" % options.save):
            log.message("") # new line
            pbar = log.pbar(len(spaths))
            import time
            for state in range(len(spaths)):
                pymol_cmd.set_frame(state + 1)
                pbar.update(state)
                time.sleep(0.1)
            pbar.finish()
            log.message("Finalizing session saving...",cont=True)  # new line
            pymol_cmd.set_frame(1)
            pymol_cmd.save(options.save, state=0)
            pymol_cmd.quit()


################################################################################


if __name__ == "__main__":
    ############################################################################
    # argument parsing
    import argparse

    description_version = '''Aqueduct library version %s
Valve driver version %s''' % (aqueduct_version_nice(), version_nice())
    description = '''Valve, Aqueduct driver'''

    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--dump-template-config", action="store_true", dest="dump_template_conf", required=False,
                        help="Dumps template config file. Suppress all other output or actions.")
    parser.add_argument("-t", action="store", dest="threads", required=False, default=None,
                        help="Limit Aqueduct calculations to given number of threads.")
    parser.add_argument("-c", action="store", dest="config_file", required=False, help="Config file filename.")
    parser.add_argument("--max-frame", action="store", dest="max_frame", required=False, help="Limit number of frames.")
    parser.add_argument("--version", action="store_true", dest="print_version", required=False,
                        help="Prints versions and exits..")

    args = parser.parse_args()
    ############################################################################
    # special option for dumping template config
    config = ValveConfig()  # config template
    if args.dump_template_conf:
        import StringIO

        config_dump = StringIO.StringIO()
        config.save_config_stream(config_dump)
        print config_dump.getvalue()
        exit(0)
    # special case of version
    if args.print_version:
        print description
        print description_version
        exit(0)

    ############################################################################
    # begin!

    valve_begin()
    valve_load_config(args.config_file, config)

    # get global options
    goptions = config.get_global_options()
    #pbar_name = goptions.pbar

    if args.threads is None:
        optimal_threads = cpu_count + 1
    else:
        optimal_threads = int(args.threads)
    log.message("Number of threads Valve is allowed to use: %d" % optimal_threads)
    if (optimal_threads > 1 and optimal_threads < 3) or (optimal_threads - 1 > cpu_count):
        log.message("Number of threads is not optimal; CPU count reported by system: %d" % cpu_count)
    # because it is used by mp.Pool it should be -1???
    if optimal_threads > 1:
        optimal_threads -= 1
        log.message("Main process would use 1 thread.")
        log.message("Concurent calculations would use %d threads." % optimal_threads)

    ############################################################################
    # STAGE 0

    reader = valve_read_trajectory(goptions.top, goptions.trj)

    if args.max_frame:
        max_frame = int(args.max_frame)
        if max_frame > reader.max_frame:
            log.warning("Desired --max-frame %d setting exceeds number of available frames (%d)." % (max_frame+1, reader.max_frame+1))
    else:
        max_frame = reader.max_frame
    log.message("Using %d of %d available frames." % (max_frame+1, reader.max_frame+1))

    # STAGE I
    result1 = valve_exec_stage(0, config, stage_I_run,
                               reader=reader,
                               max_frame=max_frame)

    # STAGE II
    result2 = valve_exec_stage(1, config, stage_II_run,
                               reader=reader,
                               max_frame=max_frame,
                               **result1)

    # STAGE III
    result3 = valve_exec_stage(2, config, stage_III_run,
                               **result2)

    # STAGE IV
    result4 = valve_exec_stage(3, config, stage_IV_run,
                               **result3)

    # STAGE V
    results = {}
    for result in (result2, result3, result4):
        results.update(result)

    result5 = valve_exec_stage(4, config, stage_V_run, no_io=True,
                               **results)

    # STAGE VI
    results = {}
    for result in (result3, result4):
        results.update(result)

    result6 = valve_exec_stage(5, config, stage_VI_run, no_io=True,
                               reader=reader,
                               **results)

    ############################################################################
    # end!

    valve_end()
