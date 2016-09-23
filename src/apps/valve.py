#!/bin/env python2
# -*- coding: utf-8 -*-

"""
This is driver for aqueduct.
"""

import logging
from aqueduct import logger, logger_name

formatter_string = '%(name)s:%(levelname)s:[%(module)s|%(funcName)s@%(lineno)d]: %(message)s'
# create and add console handler with WARNING level to the AQ logger
formatter = logging.Formatter(formatter_string)
ch = logging.StreamHandler()
ch.setFormatter(formatter)
ch.setLevel(logging.WARNING)  # default level is WARNING
logger.addHandler(ch)

import ConfigParser
import cPickle as pickle
import copy
import gzip
import multiprocessing as mp
import numpy as np
import os
import re
import operator
import shlex
import sys
from collections import namedtuple, OrderedDict
from functools import wraps
from itertools import izip_longest
from keyword import iskeyword
from scipy.spatial.distance import cdist, pdist

from multiprocessing import Manager

import MDAnalysis as mda
import roman

from aqueduct import greetings as greetings_aqueduct
from aqueduct import version as aqueduct_version
from aqueduct import version_nice as aqueduct_version_nice
from aqueduct.geom import traces
from aqueduct.geom.cluster import PerformClustering, DBSCAN, AffinityPropagation, MeanShift, KMeans, Birch
from aqueduct.geom.convexhull import is_point_within_convexhull
from aqueduct.geom.master import create_master_spath, CTypeSpathsCollection
from aqueduct.geom.smooth import WindowSmooth, MaxStepSmooth, WindowOverMaxStepSmooth, ActiveWindowSmooth, \
    ActiveWindowOverMaxStepSmooth, DistanceWindowSmooth, DistanceWindowOverMaxStepSmooth
from aqueduct.traj.dumps import TmpDumpWriterOfMDA
from aqueduct.traj.inlets import Inlets, InletTypeCodes
from aqueduct.traj.paths import GenericPaths, yield_single_paths
from aqueduct.traj.reader import ReadViaMDA
from aqueduct.traj.selections import CompactSelectionMDA, SelectionMDA
from aqueduct.utils import clui
from aqueduct.utils.helpers import range2int, Auto, what2what, lind



# TODO: Move it to separate module
cpu_count = mp.cpu_count()

# global optimal_threads
optimal_threads = None


def version():
    return 0, 9, 4


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
        options = list()
        for opt in input_options:
            if iskeyword(opt):
                logger.warning('Invalid keyword <%s> in config file skipped. Check configuration file.' % opt)
                continue
            options.append((opt, self.special_name(input_options[opt])))
        options = dict(options)
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
    def recursive_clusterization_name():
        return 'recursive_clusterization'

    @staticmethod
    def recursive_threshold_name():
        return 'recursive_threshold'

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

    def get_cluster_options(self, section_name=None):
        if section_name == None:
            section = self.cluster_name()
        else:
            section = section_name
        names = self.config.options(section)
        # options = {name: self.config.get(section, name) for name in names}
        options = dict(((name, self.config.get(section, name)) for name in names))
        # special recursive_clusterization option
        if self.recursive_clusterization_name() not in options:
            options.update({self.recursive_clusterization_name(): None})
        if self.recursive_threshold_name() not in options:
            options.update({self.recursive_threshold_name(): None})
        return self.__make_options_nt(options)

    def get_recluster_options(self):
        return self.get_cluster_options(section_name=self.recluster_name())

    def get_smooth_options(self):
        section = self.smooth_name()
        names = self.config.options(section)
        # options = {name: self.config.get(section, name) for name in names}
        options = dict(((name, self.config.get(section, name)) for name in names))
        return self.__make_options_nt(options)

    def get_default_config(self):
        # snr = 0 # stage number

        config = ConfigParser.RawConfigParser()

        def common(section):
            for setting in self.common_config_names():
                value = None
                if setting == 'execute':
                    value = 'runonce'
                elif setting == 'dump':
                    value = '%d_%s_data.dump' % (snr + 1, section)
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
        config.set(section, 'top')  # topology
        config.set(section, 'trj')  # trajectory

        # config.set(section, 'pbar', 'simple')


        ################
        snr = 0  # stage number
        # stage I
        # find traceable residues
        section = self.stage_names(snr)
        config.add_section(section)

        common(section)
        common_traj_data(section)

        ################
        snr += 1
        # stage II
        # find raw paths
        section = self.stage_names(snr)
        config.add_section(section)

        common(section)
        common_traj_data(section)

        config.set(section, 'clear_in_object_info', 'False')

        ################
        snr += 1
        # stage III
        # create separate frames
        section = self.stage_names(snr)
        config.add_section(section)

        common(section)

        config.set(section, 'auto_barber', 'False')
        config.set(section, 'discard_empty_paths', 'True')
        config.set(section, 'sort_by_id', 'True')
        config.set(section, 'apply_smoothing', 'False')
        config.set(section, 'apply_soft_smoothing', 'True')
        config.set(section, 'discard_short_paths', '1')

        ################
        snr += 1
        # stage IV
        # inlets clusterisation
        section = self.stage_names(snr)
        config.add_section(section)

        common(section)

        config.set(section, 'max_level', '5')
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
        config.set(section, self.recursive_clusterization_name(), self.cluster_name())
        config.set(section, self.recursive_threshold_name(), 'False')

        ################
        # reclusterization
        section = self.recluster_name()
        config.add_section(section)
        config.set(section, 'method', 'dbscan')
        config.set(section, self.recursive_clusterization_name(), 'False')
        config.set(section, self.recursive_threshold_name(), 'False')

        ################
        snr += 1
        # stage V
        # analysis
        section = self.stage_names(snr)
        config.add_section(section)
        common(section)
        config.remove_option(section, 'dump')
        config.set(section, 'save', '%d_%s_results.txt' % (snr + 1, section))

        config.set(section, 'dump_config', 'True')

        ################
        snr += 1
        # stage VI
        # visualize
        section = self.stage_names(snr)
        config.add_section(section)
        common(section)
        config.remove_option(section, 'dump')
        config.set(section, 'save', '%d_%s_results.py' % (snr + 1, section))

        config.set(section, 'simply_smooths', 0.05236)

        # visualize spaths, all paths in one object
        config.set(section, 'all_paths_raw', 'False')
        config.set(section, 'all_paths_smooth', 'False')
        config.set(section, 'all_paths_split', 'False')  # split by in obj out
        config.set(section, 'all_paths_raw_io', 'False')
        config.set(section, 'all_paths_smooth_io', 'False')

        # visualize spaths, separate objects
        config.set(section, 'paths_raw', 'False')
        config.set(section, 'paths_smooth', 'False')
        config.set(section, 'paths_states', 'False')
        config.set(section, 'paths_raw_io', 'False')
        config.set(section, 'paths_smooth_io', 'False')

        config.set(section, 'ctypes_raw', 'False')
        config.set(section, 'ctypes_smooth', 'False')

        # visualize clusters
        config.set(section, 'inlets_clusters', 'False')

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
        # is something missing?
        for miss in self.config.sections():
            if miss in names: continue
            output.append('[%s]' % miss)
            for option in self.config.options(miss):
                output.append('%s = %s' % (option, str(self.config.get(miss, option))))

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
        # TODO: Remove it! This is deprecated code. It, probably, never runs.
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
    return clui.gsep(sep='-', times=48)


def asep():
    return clui.gsep(sep='=', times=72)


################################################################################
# save - load helpers

def save_dump(filename, data_to_save, **kwargs):
    with clui.fbm('Saving data dump in %s file' % filename):
        with gzip.open(filename, mode='w', compresslevel=9) as f:
            # first version:
            pickle.dump({'version': version(),
                         'aqueduct_version': aqueduct_version()}, f)
            # then data to save:
            pickle.dump(data_to_save, f)
            # then other kwargs
            pickle.dump(kwargs, f)


def load_dump(filename):
    with clui.fbm('Loading data dump from %s file' % filename):
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
        logger.error('Loaded data has %s major version lower then the application.' % what)
    if current[0] < loaded[0]:
        logger.error('Loaded data has %s major version higher then the application.' % what)
    if current[0] != loaded[0]:
        logger.error('Possible problems with API compilance.')
    if current[1] > loaded[1]:
        logger.warning('Loaded data has %s minor version lower then the application.' % what)
    if current[1] < loaded[1]:
        logger.warning('Loaded data has %s minor version higher then the application.' % what)
    if current[1] != loaded[1]:
        logger.warning('Possible problems with API compilance.')


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
    assert coptions.method in ['dbscan', 'affprop', 'meanshift', 'birch',
                               'kmeans'], 'Unknown clusterization method %s.' % coptions.method

    opts = {}

    def dbscan_opts():
        if 'eps' in coptions._asdict():
            opts.update({'eps': float(coptions.eps)})
        if 'min_samples' in coptions._asdict():
            opts.update({'min_samples': int(coptions.min_samples)})
        if 'metric' in coptions._asdict():
            assert coptions.metric in ['cityblock', 'cosine', 'euclidean',
                                       'manhattan'], "Unknown metric <%s>." % coptions.metric
            opts.update({'metric': str(coptions.metric)})
        if 'algorithm' in coptions._asdict():
            assert coptions.algorithm in ['auto', 'ball_tree', 'kd_tree',
                                          'brute'], "Unknown NN algorithm <%s>." % coptions.algorithm
            opts.update({'algorithm': str(coptions.algorithm)})
        if 'leaf_size' in coptions._asdict():
            opts.update({'leaf_size': int(coptions.leaf_size)})

    def affprop_opts():
        if 'damping' in coptions._asdict():
            opts.update({'damping': float(coptions.damping)})
        if 'convergence_iter' in coptions._asdict():
            opts.update({'convergence_iter': int(coptions.convergence_iter)})
        if 'max_iter' in coptions._asdict():
            opts.update({'max_iter': int(coptions.max_iter)})
        if 'preference' in coptions._asdict():
            opts.update({'preference': float(coptions.preference)})

    def kmeans_opts():
        if 'n_clusters' in coptions._asdict():
            opts.update({'n_clusters': int(coptions.n_clusters)})
        if 'max_iter' in coptions._asdict():
            opts.update({'max_iter': int(coptions.max_iter)})
        if 'n_init' in coptions._asdict():
            opts.update({'n_init': int(coptions.n_init)})
        if 'init' in coptions._asdict():
            assert coptions.init in ['k-means++', 'random'], "Unknown initialization method <%s>." % coptions.init
            opts.update({'init': str(coptions.init)})
        if 'tol' in coptions._asdict():
            opts.update({'tol': float(coptions.tol)})

    def meanshift_opts():
        if 'cluster_all' in coptions._asdict():
            opts.update({'cluster_all': bool(coptions.cluster_all)})
        if 'bandwidth' in coptions._asdict():
            if coptions.bandwidth in (Auto, None):
                opts.update({'bandwidth': coptions.bandwidth})
            else:
                opts.update({'bandwidth': float(coptions.bandwidth)})
        if 'bin_seeding' in coptions._asdict():
            opts.update({'bin_seeding': bool(coptions.bin_seeding)})
        if 'min_bin_freq' in coptions._asdict():
            opts.update({'min_bin_freq': int(coptions.min_bin_freq)})

    def birch_opts():
        if 'threshold' in coptions._asdict():
            opts.update({'threshold': float(coptions.threshold)})
        if 'branching_factor' in coptions._asdict():
            opts.update({'branching_factor': int(coptions.branching_factor)})
        if 'n_clusters' in coptions._asdict():
            opts.update({'n_clusters': int(coptions.n_clusters)})

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
    elif coptions.method == 'birch':
        birch_opts()
        method = Birch

    return PerformClustering(method, **opts)


################################################################################

def valve_begin():
    clui.message(greetings_aqueduct())  # nice greetings
    clui.message('Aqueduct version %s' % aqueduct_version_nice())
    clui.message('Valve driver version %s' % version_nice())
    clui.message(sep())


def valve_end():
    clui.message(sep())
    clui.message('Let the Valve be always open!')
    clui.message('Goodby!')


def valve_load_config(filename, config):
    assert filename is not None, "No config file provided."
    assert os.path.isfile(filename), "Config file %s does not exist." % filename
    with clui.fbm('Load configuration file'):
        config.load_config(filename)


def valve_read_trajectory(top, traj):
    with clui.fbm('Read trajectory'):
        return TrajectoryReader(top, traj)
        # read trajectory
        # traj_list = shlex.split(traj)
        # return ReadAmberNetCDFviaMDA(top, traj_list)
        # reader = ReadDCDviaMDA(topology, trajectory)


def valve_begin_stage(stage, config):
    clui.message(sep())
    clui.message('Starting Stage %s: %s' % (roman.toRoman(stage + 1), config.stage_names(stage)))
    options = config.get_stage_options(stage)
    clui.message('Execute mode: %s' % options.execute)
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
    clui.message("Loop over frames - search of residues in object:")
    pbar = clui.pbar(max_frame)

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

    clui.message("Number of residues to trace: %d" % all_res.unique_resids_number())

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
        clui.message('Clear data on residues in object over frames.')
        clui.message('This will be recalculated on demand.')
        res_ids_in_object_over_frames = {}

    with clui.fbm("Init paths container"):
        paths = dict(((resid, GenericPaths(resid, min_pf=0, max_pf=max_frame)) for resid in
                      all_res.unique_resids(ikwid=True)))

    with reader.get() as traj_reader:

        scope = traj_reader.parse_selection(options.scope)

        with clui.fbm("Rebuild treceable residues with current trajectory"):
            all_res = rebuild_selection(all_res, traj_reader)

        clui.message("Trajectory scan:")
        pbar = clui.pbar(max_frame)

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
                is_res_in_scope = check_res_in_scope(options, scope, all_res, None)

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

    clui.message("Number of residues to trace: %d" % len(all_res))
    clui.message("Number of paths: %d" % len(paths))

    return {'all_res': all_res, 'paths': paths, 'options': options}


################################################################################

# separate_paths
def stage_III_run(config, options,
                  paths=None,
                  reader=None,
                  **kwargs):
    soptions = config.get_smooth_options()

    if options.discard_empty_paths:
        with clui.fbm("Discard residues with empty paths"):
            for key in paths.keys():
                if len(paths[key].frames) == 0:
                    paths.pop(key)

    clui.message("Create separate paths:")
    pbar = clui.pbar(len(paths))
    # yield_single_paths requires a list of paths not a dictionary
    spaths = [sp for sp, nr in yield_single_paths(paths.values(), progress=True) if pbar.update(nr + 1) is None]
    pbar.finish()

    if options.discard_short_paths > 0:
        shorter_then = int(options.discard_short_paths)
        with clui.fbm("Discard paths shorter then %d" % shorter_then):
            spaths = [sp for sp in spaths if sp.size > shorter_then]

    if options.auto_barber:
        spheres = []
        with reader.get() as traj_reader:
            clui.message("Auto Barber is looking where to cut:")
            pbar = clui.pbar(len(spaths))
            barber = traj_reader.parse_selection(options.auto_barber)
            for sp in spaths:
                if sp.has_in:
                    center = sp.coords_in[0]
                    frame = sp.path_in[0]
                    traj_reader.set_current_frame(frame)
                    radius = min(cdist(np.matrix(center), barber.atom_positions(), metric='euclidean').flatten())
                    spheres.append((center, radius))
                if sp.has_out:
                    center = sp.coords_out[-1]
                    frame = sp.path_out[-1]
                    traj_reader.set_current_frame(frame)
                    radius = min(cdist(np.matrix(center), barber.atom_positions(), metric='euclidean').flatten())
                    spheres.append((center, radius))
                pbar.update(1)
            pbar.finish()
        # remove redundant spheres
        clui.message("Removing redundant cutting places:")
        pbar = clui.pbar(len(spheres))
        some_may_be_redundant = True
        while some_may_be_redundant:
            some_may_be_redundant = False
            if len(spheres) > 1:
                for nr, sphe1 in enumerate(spheres):
                    for sphe2 in spheres[nr+1:]:
                        if sphe1[1] > sphe2[1]: continue
                        d = cdist(np.matrix(sphe1[0]), np.matrix(sphe2[0]))
                        if d + sphe1[1] < sphe2[1]:
                            spheres.pop(nr)
                            pbar.update(1)
                            some_may_be_redundant = True
                            break
                    if some_may_be_redundant:
                        break
        pbar.finish()
        clui.message("Auto Barber in action:")
        pbar = clui.pbar(len(paths))
        for p in paths.values():
            p.barber_with_spheres(spheres)
            pbar.update(1)
        pbar.finish()
        # now, it might be that some of paths are empty
        paths = {k:v for k,v in paths.iteritems() if len(v.coords) > 0}
        clui.message("Recreate separate paths:")
        pbar = clui.pbar(len(paths))
        # yield_single_paths requires a list of paths not a dictionary
        spaths = [sp for sp, nr in yield_single_paths(paths.values(), progress=True) if pbar.update(nr + 1) is None]
        pbar.finish()

        if options.discard_short_paths > 0:
            shorter_then = int(options.discard_short_paths)
            with clui.fbm("Discard (again) paths shorter then %d" % shorter_then):
                spaths = [sp for sp in spaths if sp.size > shorter_then]

    if options.sort_by_id:
        with clui.fbm("Sort separate paths by resid"):
            spaths = sorted(spaths, key=lambda sp: (sp.id.id, sp.id.nr))
    # apply smoothing?
    if options.apply_smoothing or options.apply_soft_smoothing:
        smooth = get_smooth_method(soptions)
    if options.apply_smoothing:
        clui.message('Applying hard smoothing:')
        pbar = clui.pbar(len(spaths))
        for nr, sp in enumerate(spaths):
            sp.apply_smoothing(smooth)
            pbar.update(nr + 1)
        pbar.finish()
    if options.apply_soft_smoothing:
        clui.message('Applying soft smoothing:')
        pbar = clui.pbar(len(spaths))
        for nr, sp in enumerate(spaths):
            sp.get_coords(smooth=smooth)
            pbar.update(nr + 1)
        pbar.finish()

    clui.message("Number of paths: %d" % len(paths))
    clui.message("Number of spaths: %d" % len(spaths))

    return {'paths': paths, 'spaths': spaths, 'options': options, 'soptions': soptions}


################################################################################

def get_skip_size_function(rt=None):
    if not isinstance(rt, str): return None
    assert re.compile('^[<>=]+[0-9.]+$').match(rt) is not None, "Wrong threshold definition: %s" % rt
    op = re.compile('[<>=]+')
    op = ''.join(sorted(op.findall(rt)[0]))
    vl = re.compile('[0-9.]+')
    vl = float(vl.findall(rt)[0])
    operator_dict = {'>': operator.gt,
                     '=>': operator.ge,
                     '<=': operator.le,
                     '<': operator.lt}
    assert op in operator_dict.keys(), "Unsupported operator %s in threshold %s" % (op, rt)
    return lambda size_of_cluster: operator_dict[op](vl, size_of_cluster)


def potentially_recursive_clusterization(config,
                                         clusterization_name,
                                         inlets_object,
                                         message='clusterization',
                                         deep=0,
                                         max_level=5):
    with clui.fbm("Performing %s, level %d of %d" % (message, deep, max_level), cont=False):
        logger.debug('Clustering options section: %s' % clusterization_name)
        cluster_options = config.get_cluster_options(section_name=clusterization_name)
        # TODO: Print clusterization options in a nice way!
        clustering_function = get_clustering_method(cluster_options)
        # get skip_size function according to recursive_treshold
        skip_size = get_skip_size_function(cluster_options.recursive_threshold)
        inlets_object.perform_reclustering(clustering_function, skip_outliers=True, skip_size=skip_size)
    clui.message('Number of clusters detected so far: %d' % len(inlets_object.clusters_list))
    if cluster_options.recursive_clusterization:
        deep += 1
        if deep > max_level:
            return
        return potentially_recursive_clusterization(config, cluster_options.recursive_clusterization, inlets_object,
                                                    deep=deep, max_level=max_level)


# inlets_clusterization
def stage_IV_run(config, options,
                 spaths=None,
                 **kwargs):
    coptions = config.get_cluster_options()
    rcoptions = config.get_recluster_options()
    soptions = config.get_smooth_options()

    max_level = int(options.max_level)
    assert max_level >= 0

    # new style clustering
    with clui.fbm("Create inlets"):
        inls = Inlets(spaths)
    clui.message("Number of inlets: %d" % inls.size)

    def noo():
        # returns number of outliers
        if 0 in inls.clusters_list:
            return inls.clusters.count(0)
        return 0

    if inls.size > 0:
        # ***** CLUSTERIZATION *****
        potentially_recursive_clusterization(config, config.cluster_name(), inls,
                                             message='clusterization',
                                             max_level=max_level)
        # with log.fbm("Performing clusterization"):
        #    clustering_function = get_clustering_method(coptions)
        #    inls.perform_clustering(clustering_function)
        clui.message('Number of outliers: %d' % noo())
        # ***** OUTLIERS DETECTION *****
        if options.detect_outliers:
            with clui.fbm("Detecting outliers"):
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
            clui.message('Number of clusters detected so far: %d' % len(inls.clusters_list))
            clui.message('Number of outliers: %d' % noo())
        # ***** RECLUSTERIZATION *****
        if options.recluster_outliers:
            with clui.fbm("Performing reclusterization of outliers"):
                clustering_function = get_clustering_method(rcoptions)
                # perform reclusterization
                # inls.recluster_outliers(clustering_function)
                inls.recluster_cluster(clustering_function, 0)
            clui.message('Number of clusters detected so far: %d' % len(inls.clusters_list))
            clui.message('Number of outliers: %d' % noo())
        # ***** SINGLETONS REMOVAL *****
        if options.singletons_outliers:
            with clui.fbm("Removing clusters of size %d" % int(options.singletons_outliers)):
                inls.small_clusters_to_outliers(int(options.singletons_outliers))
            clui.message('Number of clusters detected so far: %d' % len(inls.clusters_list))
            clui.message('Number of outliers: %d' % noo())

        with clui.fbm("Calculating cluster types"):
            ctypes = inls.spaths2ctypes(spaths)

        # now, there is something to do with ctypes!
        # we can create master paths!

        clui.message("Creating master paths for cluster types:")

        smooth = get_smooth_method(soptions)
        master_paths = {}
        master_paths_smooth = {}
        ctypes_generic = [ct.generic for ct in ctypes]
        ctypes_generic_list = sorted(list(set(ctypes_generic)))

        '''
        '''
        # create pool of workers - mapping function
        map_fun = map
        if optimal_threads > 1:
            pool = mp.Pool(optimal_threads)
            map_fun = pool.map


        from aqueduct.geom.master import calculate_master
        #from multiprocessing import Lock,Manager
        #lock = Manager().Lock()

        master_paths = {}
        master_paths_smooth = {}

        '''
        with clui.fbm("Calculating master and smooth paths"):
            for mpnr,mapa in enumerate(map_fun(calculate_master,[(lind(spaths, what2what(ctypes_generic, [ct])),nr,ct,None) for nr,ct in enumerate(ctypes_generic_list)] + [(lind(spaths, what2what(ctypes_generic, [ct])),nr,ct,smooth) for nr,ct in enumerate(ctypes_generic_list)])):
                if mpnr < len(ctypes_generic_list):
                    master_paths.update({ctypes_generic_list[mpnr]:mapa})
                else:
                    master_paths_smooth.update({ctypes_generic_list[mpnr-len(ctypes_generic_list)]: mapa})

        # destroy pool of workers
        if optimal_threads > 1:
            pool.close()
            pool.join()
            del pool

        '''
        pbar = clui.pbar(len(spaths)*2)

        for nr, ct in enumerate(ctypes_generic_list):
            logger.debug('CType %s (%d)' % (str(ct), nr))
            sps = lind(spaths, what2what(ctypes_generic, [ct]))
            logger.debug('CType %s (%d), number of spaths %d' % (str(ct), nr, len(sps)))
            # print len(sps),ct
            ctspc = CTypeSpathsCollection(spaths=sps,ctype=ct,pbar=pbar,threads=optimal_threads)
            master_paths.update({ct: ctspc.get_master_path(resid=nr)})
            master_paths_smooth.update({ct: ctspc.get_master_path(resid=nr, smooth=smooth)})
            del ctspc
            #master_paths.update({ct: create_master_spath(sps, resid=nr, ctype=ct, pbar=pbar)})
            #master_paths_smooth.update(
            #    {ct: create_master_spath(sps, resid=nr, ctype=ct, smooth=smooth, pbar=pbar)})
        pbar.finish()
        # TODO: issue warinig if creation of master path failed


    else:
        clui.message("No inlets found. Clusterization skipped.")
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
            clui.message(info2print)
        print >> self.filehandle, info2print

    def sep(self):
        self(asep())

    def thead(self, info2print):
        self(clui.thead(info2print))

    def under(self, info2print):
        self(clui.underline(info2print))


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
        if len(t) > 1: # traces.length_step_std requires at least 2 points
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
        clui.message('Using user provided file (%s).' % options.save)
        clui.message(sep())
        clui.message('')
    else:
        clui.message('Using standard output.')
        clui.message(sep())
        clui.message('')

    ############
    pa.sep()
    pa('Aqueduct analysis')
    pa(clui.get_str_timestamp())

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


def is_pymol_connector_session(filename):
    session_ext = re.compile('[.][pP][sS][eE]')
    if filename:
        ext = os.path.splitext(filename)[-1]
        if session_ext.match(ext):
            return True
    return False


def is_pymol_connector_script(filename):
    script = re.compile('.*[.]([pP][yY]|[.pP][yY][.][gG][zZ])$')
    if filename:
        if script.match(filename):
            return True
    return False


def stage_VI_run(config, options,
                 reader=None,
                 spaths=None,
                 inls=None,
                 ctypes=None,
                 master_paths=None,
                 master_paths_smooth=None,
                 **kwargs):
    from aqueduct.visual.pymol_connector import ConnectToPymol, SinglePathPlotter
    # from aqueduct.visual.pymol_connector import cmd as pymol_cmd
    from aqueduct.visual.helpers import ColorMapDistMap

    soptions = config.get_smooth_options()
    smooth = get_smooth_method(soptions)

    # start pymol
    with clui.fbm("Starting PyMOL"):
        pymol_connector = ConnectToPymol()
        if is_pymol_connector_script(options.save):
            pymol_connector.init_script(options.save)
        else:
            pymol_connector.init_pymol()

        if options.simply_smooths:
            spp = SinglePathPlotter(pymol_connector, linearize=float(options.simply_smooths))
        else:
            spp = SinglePathPlotter(pymol_connector, linearize=None)

    ctypes_generic = [ct.generic for ct in ctypes]
    ctypes_generic_list = sorted(list(set(ctypes_generic)))

    if options.show_molecule:
        with clui.fbm("Molecule"):
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
        with clui.fbm("Convexhull"):
            with reader.get() as traj_reader:
                options_stageII = config.get_stage_options(1)
                if options_stageII.scope_convexhull:
                    scope = traj_reader.parse_selection(options_stageII.scope)
                    frames_to_show = range2int(options.show_chull_frames)
                    for frame in frames_to_show:
                        traj_reader.set_current_frame(frame)
                        chull = scope.get_convexhull_of_atom_positions()
                        spp.convexhull(chull, state=frame + 1)

    if options.inlets_clusters:
        with clui.fbm("Clusters"):
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
        with clui.fbm("CTypes raw"):
            for nr, ct in enumerate(ctypes_generic_list):
                clui.message(str(ct), cont=True)
                sps = lind(spaths, what2what(ctypes_generic, [ct]))
                plot_spaths_traces(sps, name=str(ct) + '_raw', split=False, spp=spp)
                if ct in master_paths:
                    plot_spaths_traces([master_paths[ct]], name=str(ct) + '_raw_master', split=False, spp=spp)

    if options.ctypes_smooth:
        with clui.fbm("CTypes smooth"):
            for nr, ct in enumerate(ctypes_generic_list):
                clui.message(str(ct), cont=True)
                sps = lind(spaths, what2what(ctypes_generic, [ct]))
                plot_spaths_traces(sps, name=str(ct) + '_smooth', split=False, spp=spp, smooth=smooth)
                if ct in master_paths_smooth:
                    if master_paths_smooth[ct] is None: continue
                    plot_spaths_traces([master_paths_smooth[ct]], name=str(ct) + '_smooth_master', split=False, spp=spp,
                                       smooth=lambda anything: anything)
                if ct in master_paths:
                    if master_paths[ct] is None: continue
                    plot_spaths_traces([master_paths[ct]], name=str(ct) + '_raw_master_smooth', split=False, spp=spp,
                                       smooth=smooth)

    if options.all_paths_raw:
        with clui.fbm("All raw paths"):
            plot_spaths_traces(spaths, name='all_raw', split=options.all_paths_split, spp=spp)
    if options.all_paths_raw_io:
        with clui.fbm("All raw paths io"):
            plot_spaths_inlets(spaths, name='all_raw_paths_io', spp=spp)

    if options.all_paths_smooth:
        with clui.fbm("All smooth paths"):
            plot_spaths_traces(spaths, name='all_smooth', split=options.all_paths_split, spp=spp, smooth=smooth)
    if options.all_paths_smooth_io:
        with clui.fbm("All smooth paths io"):
            plot_spaths_inlets(spaths, name='all_smooth_paths_io', spp=spp)

    with clui.fbm("Paths as states"):
        if options.paths_raw:
            clui.message("raw", cont=True)
            plot_spaths_traces(spaths, name='raw_paths', states=options.paths_states, separate=not options.paths_states,
                               spp=spp)
        if options.paths_smooth:
            clui.message("smooth", cont=True)
            plot_spaths_traces(spaths, name='smooth_paths', states=options.paths_states,
                               separate=not options.paths_states, smooth=smooth, spp=spp)
        if options.paths_raw_io:
            clui.message("raw_io", cont=True)
            plot_spaths_inlets(spaths, name='raw_paths_io', states=options.paths_states,
                               separate=not options.paths_states, spp=spp)
        if options.paths_smooth_io:
            clui.message("smooth_io", cont=True)
            plot_spaths_inlets(spaths, name='smooth_paths_io', states=options.paths_states,
                               separate=not options.paths_states, smooth=smooth, spp=spp)

    if options.show_molecule:
        pymol_connector.orient_on('molecule')

    if is_pymol_connector_session(options.save):
        with clui.fbm("Saving session (%s)" % options.save):
            clui.message("")  # new line
            pbar = clui.pbar(len(spaths))
            import time
            for state in range(len(spaths)):
                pymol_connector.cmd.set_frame(state + 1)
                pbar.update(state)
                time.sleep(0.1)
            pbar.finish()
            clui.message("Finalizing session saving...", cont=True)  # new line
            pymol_connector.cmd.set_frame(1)
            pymol_connector.cmd.save(options.save, state=0)
            pymol_connector.cmd.quit()


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

    parser.add_argument("--debug", action="store_true", dest="debug", required=False, help="Prints debug info.")
    parser.add_argument("--debug-file", action="store", dest="debug_file", required=False, help="Debug log file.")
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
    # debug
    # at this stage logger is the AQ root logger
    if args.debug:
        logger.removeHandler(ch)  # remove old ch handlers
        ch = logging.StreamHandler()
        ch.setFormatter(formatter)
        ch.setLevel(logging.DEBUG)
        logger.addHandler(ch)
    if args.debug_file:
        formatter = logging.Formatter('%(asctime)s: ' + formatter_string)
        fh = logging.FileHandler(args.debug_file)
        fh.setFormatter(formatter)
        fh.setLevel(logging.DEBUG)
        logger.addHandler(fh)
    # finally, get valve logger
    logger = logging.getLogger(logger_name + '.valve')
    logger.info('Initialization of Valve logging done.')

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
    # pbar_name = goptions.pbar

    if args.threads is None:
        optimal_threads = cpu_count + 1
    else:
        optimal_threads = int(args.threads)
    clui.message("Number of threads Valve is allowed to use: %d" % optimal_threads)
    if (optimal_threads > 1 and optimal_threads < 3) or (optimal_threads - 1 > cpu_count):
        clui.message("Number of threads is not optimal; CPU count reported by system: %d" % cpu_count)
    # because it is used by mp.Pool it should be -1???
    if optimal_threads > 1:
        optimal_threads -= 1
        clui.message("Main process would use 1 thread.")
        clui.message("Concurent calculations would use %d threads." % optimal_threads)

    ############################################################################
    # STAGE 0

    reader = valve_read_trajectory(goptions.top, goptions.trj)

    if args.max_frame:
        max_frame = int(args.max_frame)
        if max_frame > reader.max_frame:
            logger.warning("Desired --max-frame %d setting exceeds number of available frames (%d)." % (
                max_frame + 1, reader.max_frame + 1))
    else:
        max_frame = reader.max_frame
    clui.message("Using %d of %d available frames." % (max_frame + 1, reader.max_frame + 1))

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
                               reader=reader,
                               **result2)

    # STAGE IV
    result4 = valve_exec_stage(3, config, stage_IV_run,
                               **result3)

    # STAGE V
    results = {}
    for result in (result2, result3, result4):
        results.update(result)

    result5 = valve_exec_stage(4, config, stage_V_run,
                               no_io=True,
                               **results)

    # STAGE VI
    results = {}
    for result in (result3, result4):
        results.update(result)

    result6 = valve_exec_stage(5, config, stage_VI_run,
                               no_io=True,
                               reader=reader,
                               **results)
    ############################################################################
    # end!


    ############################################################################
    # end!

    valve_end()
    logger.info('Valve calulations finished.')
