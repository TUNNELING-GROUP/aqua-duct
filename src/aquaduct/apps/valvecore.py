# -*- coding: utf-8 -*-

# This program is distributed in the hope that it will be useful,
# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2016-2018  Tomasz Magdziarz, Alicja Płuciennik, Michał Stolarczyk <info@aquaduct.pl>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


#import ipdb as pdb


import ConfigParser
import copy
import multiprocessing as mp
import numpy as np
import operator
import os
import re
import shlex
import sys
from collections import namedtuple, OrderedDict
from functools import wraps, partial
from itertools import izip_longest, izip, imap
from keyword import iskeyword

import array

from scipy.spatial.distance import cdist
from scipy.stats import ttest_ind

from aquaduct.utils.clui import roman

from aquaduct import greetings as greetings_aquaduct
from aquaduct import logger
from aquaduct import version_nice as aquaduct_version_nice
from aquaduct.apps.data import get_vda_reader, GCS, CRIC, save_cric
from aquaduct.geom import traces
from aquaduct.geom.cluster import AVAILABLE_METHODS as available_clusterization_methods
from aquaduct.geom.cluster import PerformClustering, DBSCAN, AffinityPropagation, MeanShift, KMeans, Birch, \
    BarberCluster, get_required_params
from aquaduct.geom.master import CTypeSpathsCollection
from aquaduct.geom.smooth import WindowSmooth, MaxStepSmooth, WindowOverMaxStepSmooth, ActiveWindowSmooth, \
    ActiveWindowOverMaxStepSmooth, DistanceWindowSmooth, DistanceWindowOverMaxStepSmooth, SavgolSmooth
from aquaduct.traj.barber import WhereToCut
from aquaduct.traj.dumps import TmpDumpWriterOfMDA
from aquaduct.traj.inlets import InletClusterGenericType
from aquaduct.traj.inlets import Inlets, InletTypeCodes
from aquaduct.traj.paths import GenericPaths, yield_single_paths, PassingPath, SinglePath
from aquaduct.traj.paths import union_full, yield_generic_paths
# from aquaduct.traj.reader import ReadViaMDA
# from aquaduct.traj.selections import CompactSelectionMDA
from aquaduct.utils import clui
from aquaduct.utils.helpers import range2int, Auto, what2what, lind, is_number, robust_and, robust_or
from aquaduct.utils.multip import optimal_threads
from aquaduct.traj.sandwich import ResidueSelection, Reader
from aquaduct.utils.helpers import SmartRange, iterate_or_die

__mail__ = 'info@aquaduct.pl'
__version__ = aquaduct_version_nice()


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


# TODO: add parser for initial checking of configuration file
class ValveConfig(object, ConfigSpecialNames):
    def __init__(self):
        self.config = self.get_default_config()
        self.config_filename = ''

    def __make_options_nt(self, input_options):
        # options = {opt: self.special_name(input_options[opt]) for opt in input_options}
        options = list()
        for opt in input_options:
            if iskeyword(opt):
                logger.warning('Invalid keyword <%s> in config file skipped. Check configuration file.' % opt)
                continue
            if opt.replace('_', '').isalnum():
                options.append((opt, self.special_name(input_options[opt])))
            else:
                logger.debug('Keyword <%s> in config file skipped.' % opt)
        options = OrderedDict(options)
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
        return 'scope scope_convexhull scope_everyframe object'.split()

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
        options = OrderedDict(((name, None) for name in self.common_traj_data_config_names()))
        for nr in range(stage + 1)[::-1]:
            section = self.stage_names(nr)
            for name in self.common_traj_data_config_names():
                if self.config.has_option(section, name):
                    value = self.special_name(self.config.get(section, name))
                    if (value is not None) and (options[name] is None):
                        options.update({name: value})
        return self.__make_options_nt(options)

    def get_global_options(self):
        section = self.global_name()
        names = self.config.options(section)
        # options = {name: self.config.get(section, name) for name in names}
        options = OrderedDict(((name, self.config.get(section, name)) for name in names))
        return self.__make_options_nt(options)

    def get_stage_options(self, stage):
        assert isinstance(stage, int)
        stage_name = self.stage_names(stage)
        names = self.config.options(stage_name)
        # options = {name: self.config.get(stage_name, name) for name in names}
        options = OrderedDict(((name, self.config.get(stage_name, name)) for name in names))
        if stage in [0, 1]:
            options.update(self.get_common_traj_data(stage)._asdict())
        return self.__make_options_nt(options)

    def get_cluster_options(self, section_name=None):
        if section_name is None:
            section = self.cluster_name()
        else:
            section = section_name
        names = self.config.options(section)
        # options = {name: self.config.get(section, name) for name in names}
        options = OrderedDict(((name, self.config.get(section, name)) for name in names))
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
        options = OrderedDict(((name, self.config.get(section, name)) for name in names))
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

        def common_traj_data(section, commented=False):
            for setting in self.common_traj_data_config_names():
                if commented:
                    setting = '#' + setting
                config.set(section, setting, 'None')

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
        config.set(section, 'scope_everyframe', 'False')

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

        config.set(section, 'allow_passing_paths', 'False')
        config.set(section, 'auto_barber_tovdw', 'True')
        config.set(section, 'auto_barber_maxcut', '2.8')
        config.set(section, 'auto_barber_mincut', 'None')
        config.set(section, 'auto_barber_maxcut_level', 'True')
        config.set(section, 'auto_barber_mincut_level', 'True')
        config.set(section, 'auto_barber', 'False')
        config.set(section, 'discard_empty_paths', 'True')
        config.set(section, 'sort_by_id', 'True')
        config.set(section, 'apply_smoothing', 'False')
        config.set(section, 'apply_soft_smoothing', 'True')

        config.set(section, 'discard_short_paths', '20')
        config.set(section, 'discard_short_object', '2.0')
        config.set(section, 'discard_short_logic', 'or')

        ################
        snr += 1
        # stage IV
        # inlets clusterisation
        section = self.stage_names(snr)
        config.add_section(section)

        common(section)

        config.set(section, 'max_level', '0')
        config.set(section, 'recluster_outliers', 'False')
        config.set(section, 'detect_outliers', 'False')
        config.set(section, 'singletons_outliers', 'False')
        config.set(section, 'create_master_paths', 'False')
        config.set(section, 'exclude_passing_in_clusterization', 'True')
        config.set(section, 'add_passing_to_clusters', 'None')

        ################
        # smooth
        section = self.smooth_name()
        config.add_section(section)
        config.set(section, 'method', 'window')

        ################
        # clusterization
        section = self.cluster_name()
        config.add_section(section)
        config.set(section, 'method', 'barber')
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

        config.set(section, 'calculate_scope_object_size', 'False')
        config.set(section, 'scope_chull', 'None')
        config.set(section, 'object_chull', 'None')

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

        config.set(section, 'simply_smooths', 'RecursiveVector')
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
        config.set(section, 'show_scope_chull', 'None')
        config.set(section, 'show_scope_chull_frames', '0')
        config.set(section, 'show_object_chull', 'None')
        config.set(section, 'show_object_chull_frames', '0')

        return config

    def load_config(self, filename):
        self.config.read(filename)
        self.config_filename = filename

    def save_config_stream(self, fs):
        self.config.write(fs)

    def save_config(self, filename):
        with open(filename, 'w') as fs:
            self.save_config_stream(fs)

    def get_general_comment(self, section):
        if section == self.global_name():
            return ['Global settings.']
        if section == self.cluster_name():
            out = ['Possible clusterization methods:']
            # get default clusterization method
            defmet = self.config.get(section, 'method')
            for method in available_clusterization_methods:
                if method == defmet: continue
                out.append('method = %s' % method)
                params = get_required_params(method)
                if params:
                    sss = ''
                    if len(params) > 1:
                        sss = 's'
                    out.append('Required parameter%s for %s method:' % (sss, method))
                    for param in params:
                        out.append('%s = None' % param)
                    # out.append('')
            return out

    def dump_config(self, dump_template=False):
        skip_list = '''simply_smooths'''.split()  # these options are very optional!
        concise = True
        if dump_template:
            concise = False

        output = []

        options = [self.get_global_options()] + \
                  [self.get_stage_options(stage) for stage in range(6)] + \
                  [self.get_cluster_options(), self.get_recluster_options(),
                   self.get_smooth_options()]
        names = [self.global_name()] + \
                [self.stage_names(stage) for stage in range(6)] + \
                [self.cluster_name(), self.recluster_name(),
                 self.smooth_name()]

        def value2str(val):
            if val is Auto:
                val = "Auto"
            else:
                val = str(val)
            return val

        # this loop ensures it is returned in particular order
        for opts, name in zip(options, names):
            output.append('[%s]' % name)  # section name
            # general comments
            if dump_template:
                comment = self.get_general_comment(name)
                if comment:
                    output.extend(['# %s' % line for line in comment])
            for key, value in opts._asdict().iteritems():  # loop over options
                if key in skip_list: continue
                # comment scope etc. in stage II if dump_template
                if dump_template and name == self.stage_names(1):
                    if key in self.common_traj_data_config_names():
                        key = '#' + key
                output.append('%s = %s' % (key, value2str(value)))
            if not concise: output.append('')

        # is something missing? another loop over all additional sections
        for miss in self.config.sections():
            if miss in names: continue  # skip if it was already dumped in the loop above
            output.append('[%s]' % miss)
            for key in self.config.options(miss):  # loop over options in section miss
                if key in skip_list: continue
                value = self.config.get(miss, key)
                output.append('%s = %s' % (key, value2str(value)))
            if not concise: output.append('')
        while not len(output[-1].strip()):
            output.pop()
        return output


################################################################################

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

def get_smooth_method(soptions):
    assert soptions.method in ['window', 'mss', 'window_mss', 'awin',
                               'awin_mss', 'dwin', 'dwin_mss',
                               'savgol'], 'Unknown smoothing method %s.' % soptions.method

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

    def savgol_opts():
        window_length = 5  # TODO: magic constant (default value)
        if 'window' in soptions._asdict():
            window_length = int(float(soptions.window))
            assert window_length % 2 == 1, 'Window in Savgol method must be positive odd number, %d given instead.' % window_length
            opts.update({'window_length': window_length})
        polyorder = 2  # TODO: magic constant (default value)
        if 'polyorder' in soptions._asdict():
            polyorder = int(float(soptions.polyorder))
            assert polyorder > 0, 'Polynomial order should be greater then 0, %d given instead.' % polyorder
            opts.update({'polyorder': polyorder})
        assert polyorder < window_length, 'Polynomial order (%d) should be less then window (%d).' % (
            polyorder, window_length)
        '''
        if 'deriv' in soptions._asdict():
            deriv = int(float(soptions.deriv))
            assert deriv >= 0 and deriv <= polyorder, 'Order of derrivative should be integer greater or equal 0 and less then or equal to Polynomial order (%d), %d given instead.' % (polyorder,deriv)
            opts.update({'deriv': deriv})
        if 'delta' in soptions._asdict():
            delta = float(soptions.delta)
            assert delta >= 0, 'Delta should be greater or equal 0, %f given instead.' % delta
            if 'deriv' in soptions._asdict():
                if deriv == 0:
                    logger.warning('Delta %f make no sense if deriv is 0.' % delta)
            opts.update({'delta': delta})
        mode = 'interp' # TODO: magic constant (default value)
        if 'mode' in soptions._asdict():
            mode = str(soptions.mode)
            assert mode in ['mirror', 'constant', 'nearest', 'wrap', 'interp'], 'Unknown mode %s.' % mode
            opts.update({'mode': mode})
        if 'cval' in soptions._asdict():
            cval = float(soptions.cval)
            if mode != 'constant':
                logger.warning('Cval make no sense if mode is %s.' % mode)
            opts.update({'cval': cval})
        '''

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
    elif soptions.method == 'savgol':
        savgol_opts()
        smooth = SavgolSmooth(**opts)

    return smooth


def get_auto_barber_options(abo):
    opts = {}
    if 'auto_barber' in abo._asdict():
        opts.update({'selection': str(abo.auto_barber)})
    if 'auto_barber_maxcut' in abo._asdict():
        if is_number(abo.auto_barber_maxcut):
            opts.update({'maxcut': float(abo.auto_barber_maxcut)})
        else:
            opts.update({'maxcut': None})
    if 'auto_barber_mincut' in abo._asdict():
        if is_number(abo.auto_barber_mincut):
            opts.update({'mincut': float(abo.auto_barber_mincut)})
        else:
            opts.update({'mincut': None})
    if 'auto_barber_mincut_level' in abo._asdict():
        opts.update({'mincut_level': bool(abo.auto_barber_mincut_level)})
    if 'auto_barber_maxcut_level' in abo._asdict():
        opts.update({'maxcut_level': bool(abo.auto_barber_maxcut_level)})
    if 'auto_barber_tovdw' in abo._asdict():
        opts.update({'tovdw': bool(abo.auto_barber_tovdw)})
    return opts


def get_clustering_method(coptions, config):
    assert coptions.method in available_clusterization_methods, 'Unknown clusterization method %s.' % coptions.method

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

    def barber_opts():
        abo = config.get_stage_options(2)
        opts.update(get_auto_barber_options(abo))
        opts.update(get_auto_barber_options(coptions))

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
    elif coptions.method == 'barber':
        barber_opts()
        method = BarberCluster

    return PerformClustering(method, **opts)


def get_linearize_method(loption):
    if loption:
        assert isinstance(loption, (str, unicode)), "Wrong Linearize method definition: %r" % loption
        possible_formats = [
            re.compile('^(recursive|oneway|hobbit)(triangle|vector)[(][+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?[)]$'),
            re.compile('^(recursive|oneway|hobbit)(triangle|vector)[(][)]$'),
            re.compile('^(recursive|oneway|hobbit)(triangle|vector)$')]
        assert True in [pf.match(loption.lower()) is not None for pf in
                        possible_formats], "Wrong Linearize method definition: %s" % loption
        # http://stackoverflow.com/questions/12929308/python-regular-expression-that-matches-floating-point-numbers#12929311
        way = [w for w in ['recursive', 'oneway', 'hobbit'] if w in loption.lower()][0]
        crit = [c for c in ['triangle', 'vector'] if c in loption.lower()][0]
        threshold = re.compile('[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?').findall(loption)
        if len(threshold):
            threshold = float(threshold[0][0])
        else:
            threshold = None
        # get method
        if way == 'recursive':
            if crit == 'triangle':
                met = traces.LinearizeRecursiveTriangle
            else:
                met = traces.LinearizeRecursiveVector
        elif way == 'oneway':
            if crit == 'triangle':
                met = traces.LinearizeOneWayTriangle
            else:
                met = traces.LinearizeOneWayVector
        elif way == 'hobbit':
            if crit == 'triangle':
                met = traces.LinearizeHobbitTriangle
            else:
                met = traces.LinearizeHobbitVector
        if threshold is None:
            return met()
        return met(threshold)


################################################################################

def valve_begin():
    clui.message(greetings_aquaduct())  # nice greetings
    clui.message('Aqua-Duct version %s' % aquaduct_version_nice())
    # clui.message('Valve driver version %s' % version_nice())
    clui.message(sep())


def valve_end():
    clui.message(sep())
    clui.message('Let the Valve be always open!')
    clui.message('Goodbye!')


def valve_load_config(filename, config):
    assert filename is not None, "No config file provided."
    assert os.path.isfile(filename), "Config file %s does not exist." % filename
    with clui.fbm('Load configuration file'):
        config.load_config(filename)


'''
def valve_read_trajectory(top, traj, frames_window=None):
    with clui.fbm('Read trajectory'):
        return TrajectoryReader(top, traj, frames_window=frames_window)
        # read trajectory
        # traj_list = shlex.split(traj)
        # return ReadAmberNetCDFviaMDA(top, traj_list)
        # reader = ReadDCDviaMDA(topology, trajectory)
'''


def valve_begin_stage(stage, config):
    clui.message(sep())
    clui.message('Starting Stage %s: %s' % (roman.toRoman(stage + 1), config.stage_names(stage)))
    options = config.get_stage_options(stage)
    clui.message('Execute mode: %s' % options.execute)
    return options


def valve_exec_stage(stage, config, stage_run, no_io=False, run_status=None,
                     **kwargs):
    with clui.tictoc('Stage %s (%s)' % (roman.toRoman(stage + 1), config.stage_names(stage))):

        # This function runs stages in a smart way, checks execution logic, and loads/saves dumps if required.
        options = valve_begin_stage(stage, config)

        run_status.update({stage: False})

        # TODO: Consider to create traj_reader object here instead of doing it in stage_run or in load...
        # execute?
        can_be_loaded = False
        if (not no_io) and options.dump:
            if os.path.isfile(options.dump) or os.path.islink(options.dump):
                can_be_loaded = True
        # has to be run?
        if options.execute in ['runonce'] and can_be_loaded and stage > 0:
            if run_status[stage - 1]:
                can_be_loaded = False

        if options.execute in ['run'] or (options.execute in ['runonce'] and not can_be_loaded):
            result = stage_run(config, options, **kwargs)
            save_cric()
            run_status.update({stage: True})
            if not no_io:
                ###########
                # S A V E #
                ###########
                with clui.fbm('Saving data dump in %s file' % options.dump):
                    vda = get_vda_reader(options.dump)
                    vda(mode='w', data_file_name=options.dump).dump(**result)
                # save_stage_dump(options.dump, **result)
        elif options.execute in ['skip'] or (options.execute in ['runonce'] and can_be_loaded):
            if not no_io:
                ###########
                # L O A D #
                ###########
                if options.dump:
                    with clui.fbm('Loading data dump from %s file' % options.dump):
                        vda = get_vda_reader(options.dump)
                        result = vda(mode='r', data_file_name=options.dump).load()
                        # result = load_stage_dump(options.dump, reader=reader)
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
                **kwargs):
    # create pool of workers - mapping function
    map_fun = map
    if optimal_threads.threads_count > 1:
        pool = mp.Pool(optimal_threads.threads_count)
        map_fun = pool.map

    clui.message("Loop over frames - search of residues in object:")
    pbar = clui.pbar(Reader.number_of_frames())

    # create some containers
    number_frame_rid_in_object = []
    # res_ids_in_scope_over_frames = {}  # not used
    all_res = None

    # scope will be used to derrive center of system
    center_of_system = np.array([0., 0., 0.])

    # loop over possible layers of sandwich
    for number, traj_reader in Reader.iterate(number=True):

        # scope is evaluated only once before the loop over frames starts
        if not options.scope_everyframe:
            scope = traj_reader.parse_selection(options.scope)

        # all_res_this_layer = [] # list of all new res in this layer found in an order of apparance

        frame_rid_in_object = []

        # the loop over frames
        for frame in traj_reader.iterate_over_frames():
            if options.scope_everyframe:
                scope = traj_reader.parse_selection(options.scope)
            # center of system
            center_of_system += scope.center_of_mass()
            # current res selection
            res = traj_reader.parse_selection(options.object).residues()
            # find matching residues, ie those which are in the scope:
            res_new = scope.containing_residues(res, convex_hull=options.scope_convexhull, map_fun=map_fun)
            res_new.uniquify()  # here is a list of residues in this layer that are in the object and in the scope
            # adds them to all_res
            if all_res:
                all_res.add(res_new)
                all_res.uniquify()
            else:
                all_res = res_new
            # remeber ids of res in object in current frame
            if res_new is not None:
                frame_rid_in_object.append([rid[-1] for rid in res_new.ids()])
            else:
                frame_rid_in_object.append([])
            pbar.next()

        number_frame_rid_in_object.append(frame_rid_in_object)

    # destroy pool of workers
    if optimal_threads.threads_count > 1:
        pool.close()
        pool.join()
        del pool

    # res_ids_in_object_frames_list = [res_ids_in_object_frames_list[number][rid] for number,rid in all_res.ids()]

    pbar.finish()
    center_of_system /= (Reader.number_of_frames())
    logger.info('Center of system is %0.2f, %0.2f, %0.2f' % tuple(center_of_system))

    if all_res is None:
        raise ValueError("No traceable residues was found.")

    clui.message("Number of residues to trace: %d" % all_res.len())

    # 'res_ids_in_object_over_frames': IdsOverIds.dict2arrays(res_ids_in_object_over_frames),
    return {'all_res': all_res,
            # 'res_ids_in_object_over_frames': res_ids_in_object_over_frames,
            'number_frame_rid_in_object': number_frame_rid_in_object,
            'center_of_system': center_of_system,
            'options': options._asdict()}


################################################################################

# raw_paths
def stage_II_run(config, options,
                 all_res=None,
                 number_frame_rid_in_object=None,
                 # res_ids_in_object_over_frames=None,
                 **kwargs):
    ####################################################################################################################
    # FIXME: temporary solution, remove it later
    if 'res_ids_in_object_over_frames' in kwargs:
        res_ids_in_object_over_frames = kwargs['res_ids_in_object_over_frames']
    else:
        res_ids_in_object_over_frames = None
    if number_frame_rid_in_object is None and res_ids_in_object_over_frames is not None:
        with clui.fbm("Get number_frame_rid_in_object if possible"):
            number_frame_rid_in_object = []
            for number in xrange(max(res_ids_in_object_over_frames.keys())):
                frame_rid_in_object = []
                for layer in xrange(max(res_ids_in_object_over_frames[number].keys())):
                    frame_rid_in_object.append(res_ids_in_object_over_frames[number][layer])
                number_frame_rid_in_object.append(frame_rid_in_object)
    ####################################################################################################################
    # create pool of workers - mapping function
    map_fun = map
    if optimal_threads.threads_count > 1:
        pool = mp.Pool(optimal_threads.threads_count)
        chunk_size = all_res.len() / Reader.number_of_layers() / (optimal_threads.threads_count ** 1) - 1
        if chunk_size <= 0:
            chunk_size = 1
        logger.debug("Chunk size %d.",chunk_size)
        map_fun = partial(pool.imap, chunksize=chunk_size)
        #map_fun = pool.map
    # clear in object info if required
    if options.clear_in_object_info:
        clui.message('Clear data on residues in object over frames.')
        clui.message('This will be recalculated on demand.')
        number_frame_rid_in_object = None

    is_number_frame_rid_in_object = bool(number_frame_rid_in_object)

    '''
    with clui.fbm("Init paths container"):
        number_of_frames = Reader.window.len() - 1
        paths = dict(
            ((resid, GenericPaths(resid, name_of_res=resname, single_res_selection=sressel,
                                  min_pf=0, max_pf=number_of_frames))
             for resid, resname, sressel in
             zip(all_res.ids(), all_res.names(), all_res.single_residues())))
    '''
    number_of_frames = Reader.number_of_frames(onelayer=True) - 1
    clui.message("Trajectory scan:")
    pbar = clui.pbar(Reader.number_of_frames())
    # loop over possible layers of sandwich
    paths = []
    for frame_rid_in_object, (number, traj_reader) in izip(
            iterate_or_die(number_frame_rid_in_object, times=Reader.number_of_layers()), Reader.iterate(number=True)):

        # scope is evaluated only once before loop over frames so it cannot be frame dependent
        if not options.scope_everyframe:
            scope = traj_reader.parse_selection(options.scope)
            logger.debug("Scope definition evaluated only once for given layer")
        else:
            logger.debug("Scope definition evaluated in every frame, this might be very slow.")

        # speed up!
        all_res_this_layer = all_res.layer(number)
        all_res_this_ids = list(all_res_this_layer.ids())

        paths_this_layer = (GenericPaths(resid,
                                         name_of_res=resname,
                                         single_res_selection=sressel,
                                         min_pf=0, max_pf=number_of_frames)
                            for resid, resname, sressel in izip(all_res_this_ids,
                                                                all_res_this_layer.names(),
                                                                all_res_this_layer.single_residues()))

        # big container for 012 path data
        number_frame_object_scope = np.zeros((Reader.number_of_frames(onelayer=True), all_res_this_layer.len()), dtype=np.int8)
        # the loop over frames, use izip otherwise iteration over frames does not work
        for rid_in_object, frame in izip(
                iterate_or_die(frame_rid_in_object, times=Reader.number_of_frames(onelayer=True)),
                traj_reader.iterate_over_frames()):

            # do we have object data?
            if not is_number_frame_rid_in_object:
                rid_in_object = [rid[-1] for rid in traj_reader.parse_selection(options.object).residues().ids()]
            # assert rid_in_object is not None

            is_res_in_object = (rid[-1] in rid_in_object for rid in all_res_this_ids)

            if options.scope_everyframe:
                scope = traj_reader.parse_selection(options.scope)
            # check if all_res are in the scope, reuse res_ids_in_object_over_frames
            is_res_in_scope = scope.contains_residues(all_res_this_layer, convex_hull=options.scope_convexhull,
                                                      map_fun=map_fun,
                                                      known_true=None)  # known_true could be rid_in_object

            number_frame_object_scope[frame,:] = np.array(map(sum,izip(is_res_in_object,is_res_in_scope)),dtype=np.int8)
            pbar.next()

        #number_frame_object_scope = np.array(number_frame_object_scope,dtype=np.int8).T
        # another loop over this columns
        for pat,nfos in izip(paths_this_layer,number_frame_object_scope.T):
            pat.add_012(nfos)
            paths.append(pat)
        del number_frame_object_scope

        #paths.extend(paths_this_layer)

    # destroy pool of workers
    if optimal_threads.threads_count > 1:
        pool.close()
        pool.join()
        del pool

    pbar.finish()

    clui.message("Number of paths: %d" % len(paths))

    return {'all_res': all_res, 'paths': paths, 'options': options._asdict()}


################################################################################

class ABSphere(namedtuple('ABSphere', 'center radius')):
    pass


# separate_paths
def stage_III_run(config, options,
                  paths=None,
                  **kwargs):
    soptions = config.get_smooth_options()

    if options.allow_passing_paths:
        logger.warning("Passing paths is a highly experimental feature. Please, analyze results with care.")

    if options.discard_empty_paths:
        with clui.fbm("Discard residues with empty paths"):
            paths = [pat for pat in paths if len(pat.frames) > 0]

    clui.message("Create separate paths:")
    pbar = clui.pbar(len(paths))
    # yield_single_paths requires a list of paths not a dictionary
    spaths = [sp for sp, nr in yield_single_paths(paths,
                                                  progress=True,
                                                  passing=options.allow_passing_paths) if pbar.update(nr + 1) is None]
    pbar.finish()
    clui.message("Created %d separate paths out of %d raw paths" %
                 (len(spaths),len(paths)))
    pbar = clui.pbar(len(spaths),"Removing unused parts of paths:")
    paths = yield_generic_paths(spaths,progress=pbar)
    pbar.finish()

    if options.discard_short_paths or options.discard_short_object:
        if is_number(options.discard_short_paths):
            short_paths = int(options.discard_short_paths)
        else:
            short_paths = None
        if is_number(options.discard_short_object):
            short_object = float(options.discard_short_object)
        else:
            short_object = None

        if options.discard_short_logic == 'and':
            short_logic = robust_or  # for and use or, this is intentional
            short_logic_name = "AND"
        else:
            short_logic = robust_and
            short_logic_name = "OR"
            if options.discard_short_logic != 'or':
                logger.warning("Invalid discard_short_logic '%s', using %s by default." % (
                options.discard_short_logic, short_logic_name))
        # make message
        if short_paths is not None and short_object is not None:
            discard_message = "Discard paths shorter than %d %s object shorter than %0.2f" % (
            short_paths, short_logic_name, short_object)
        elif short_paths is None and short_object is not None:
            discard_message = "Discard paths object shorter than %0.2f" % short_object
        elif short_paths is not None and short_object is None:
            discard_message = "Discard paths shorter than %d" % short_paths
        if short_paths is not None or short_object is not None:
            with clui.fbm(discard_message):
                spaths_nr = len(spaths)
                # TODO: if not short object is used there is no sense in calling object_len as it is very expensive
                if short_object is not None:
                    spaths = [sp for sp in spaths if short_logic(sp.size > short_paths, sp.object_len > short_object)]
                else:
                    spaths = [sp for sp in spaths if sp.size > short_paths]
                spaths_nr_new = len(spaths)
            if spaths_nr == spaths_nr_new:
                clui.message("No paths were discarded.")
            else:
                clui.message("%d paths were discarded." % (spaths_nr - spaths_nr_new))
        else:
            clui.message("No paths were discarded - no values were set.")

    if options.auto_barber:
        wtc = WhereToCut(spaths=spaths,
                         selection=options.auto_barber,
                         mincut=options.auto_barber_mincut,
                         mincut_level=options.auto_barber_mincut_level,
                         maxcut=options.auto_barber_maxcut,
                         maxcut_level=options.auto_barber_maxcut_level,
                         tovdw=options.auto_barber_tovdw)
        # cut thyself!
        wtc.cut_thyself()

        clui.message("Auto Barber in action:")
        pbar = clui.pbar(len(paths))
        for p in paths:
            p.barber_with_spheres(wtc.spheres)
            pbar.next()
        pbar.finish()
        # now, it might be that some of paths are empty
        # paths = {k: v for k, v in paths.iteritems() if len(v.coords) > 0}
        # paths = dict((k, v) for k, v in paths.iteritems() if
        #             len(v.frames) > 0)  # more universal as dict comprehension may not work in <2.7
        paths = [pat for pat in paths if len(pat.frames) > 0]

        clui.message("Recreate separate paths:")
        pbar = clui.pbar(len(paths))
        # yield_single_paths requires a list of paths not a dictionary
        spaths = [sp for sp, nr in yield_single_paths(paths, progress=True,
                                                      passing=options.allow_passing_paths) if
                  pbar.update(nr + 1) is None]
        pbar.finish()

        if options.discard_short_paths or options.discard_short_object:
            if short_paths is not None or short_object is not None:
                with clui.fbm(discard_message):
                    spaths_nr = len(spaths)
                    spaths = [sp for sp in spaths if short_logic(sp.size > short_paths, sp.object_len > short_object)]
                    spaths_nr_new = len(spaths)
                if spaths_nr == spaths_nr_new:
                    clui.message("No paths were discarded.")
                else:
                    clui.message("%d paths were discarded." % (spaths_nr - spaths_nr_new))
            else:
                clui.message("No paths were discarded - no values were set.")

    if options.sort_by_id:
        with clui.fbm("Sort separate paths by resid"):
            spaths = sorted(spaths, key=lambda sp: (sp.id.id, sp.id.nr))
    # apply smoothing?
    # it is no longer necessary
    if options.apply_smoothing:
        logger.warning("Hard smoothing is not available in the current version but may be available in the future. Stay tuned!")
    if options.apply_soft_smoothing:
        logger.warning("Soft smoothing option is not available any more. Soft smoothing is enabled by default if --cache-dir or --cache-mem options are used.")
    clui.message("Number of paths: %d" % len(paths))
    clui.message("Number of spaths: %d" % len(spaths))

    return {'paths': paths, 'spaths': spaths, 'options': options._asdict(), 'soptions': soptions._asdict()}


################################################################################

def get_allow_size_function(rt=None):
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
    return lambda size_of_cluster: operator_dict[op](size_of_cluster, vl)


class SkipSizeFunction(object):

    def __init__(self, ths_def):

        self.thresholds = []
        if isinstance(ths_def, (str, unicode)):
            for thd in ths_def.split():
                self.thresholds.append(get_allow_size_function(thd))

    def __call__(self, size_of_cluster):
        for thd in self.thresholds:
            if not thd(size_of_cluster):
                return True
        return False


def potentially_recursive_clusterization(config=None,
                                         clusterization_name=None,
                                         inlets_object=None,
                                         spaths=None,
                                         message='clusterization',
                                         deep=0,
                                         max_level=5):
    with clui.fbm("Performing %s, level %d of %d" % (message, deep, max_level), cont=False):
        logger.debug('Clustering options section: %s' % clusterization_name)
        cluster_options = config.get_cluster_options(section_name=clusterization_name)
        clui.message('Clustering options:')
        for k, v in cluster_options._asdict().iteritems():
            clui.message("%s = %s" % (str(k), str(v)))
        # TODO: Print clusterization options in a nice way!
        clustering_function = get_clustering_method(cluster_options, config)
        # special case of barber!!!
        if cluster_options.method == 'barber':
            logger.debug('Getting inltets refs...')
            inlets_refs = inlets_object.get_inlets_references()
            logger.debug('Starting wtc...')
            wtc = WhereToCut(spaths=(sp for sp in spaths if sp.id in inlets_refs),
                             expected_nr_of_spaths=len(inlets_refs),
                             forceempty=True,
                             **clustering_function.method_kwargs)
            # clouds = wtc.cloud_groups(progress=True)
            logger.debug('Getting spheres...')
            inlets_object.add_spheres(wtc.spheres)
        logger.debug('Proceed with clusterization, skip size...')
        # get skip_size function according to recursive_treshold
        skip_size = SkipSizeFunction(cluster_options.recursive_threshold)
        logger.debug('Proceed with clusterization, call method...')
        inlets_object.perform_reclustering(clustering_function, skip_outliers=True, skip_size=skip_size)
    clui.message('Number of clusters detected so far: %d' % len(inlets_object.clusters_list))

    if cluster_options.recursive_clusterization:
        deep += 1
        if deep > max_level:
            return
        return potentially_recursive_clusterization(config=config,
                                                    clusterization_name=cluster_options.recursive_clusterization,
                                                    inlets_object=inlets_object,
                                                    spaths=spaths,
                                                    deep=deep,
                                                    max_level=max_level)


# inlets_clusterization
def stage_IV_run(config, options,
                 spaths=None,
                 center_of_system=None,
                 **kwargs):
    coptions = config.get_cluster_options()
    rcoptions = config.get_recluster_options()
    soptions = config.get_smooth_options()

    max_level = int(options.max_level)
    assert max_level >= 0

    # new style clustering
    #with clui.fbm("Create inlets"):
    # here we can check center of system
    pbar = clui.SimpleProgressBar(maxval=len(spaths),mess="Create inlets")
    inls = Inlets(spaths, center_of_system=center_of_system, passing=not options.exclude_passing_in_clusterization, pbar=pbar)
    pbar.finish()
    clui.message("Number of inlets: %d" % inls.size)

    def noo():
        # returns number of outliers
        if 0 in inls.clusters_list:
            return inls.clusters.count(0)
        return 0

    if inls.size > 0:
        # ***** CLUSTERIZATION *****
        with clui.fbm("Performing clusterization", cont=False):
            potentially_recursive_clusterization(config=config,
                                                 clusterization_name=config.cluster_name(),
                                                 inlets_object=inls,
                                                 spaths=spaths,
                                                 message='clusterization',
                                                 max_level=max_level)
        # with log.fbm("Performing clusterization"):
        #    clustering_function = get_clustering_method(coptions)
        #    inls.perform_clustering(clustering_function)
        clui.message('Number of outliers: %d' % noo())
        # ***** OUTLIERS DETECTION *****
        if options.detect_outliers:
            with clui.fbm("Detecting outliers", cont=False):
                if options.detect_outliers is not Auto:
                    threshold = float(options.detect_outliers)
                clusters = list(inls.clusters)
                no_out_detected = 0
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
                            no_out_detected += 1
                clui.message('Detected %d outliers.' % no_out_detected)
                inls.add_outliers_annotations(clusters)
            clui.message('Number of clusters detected so far: %d' % len(inls.clusters_list))
            clui.message('Number of outliers: %d' % noo())
        # ***** RECLUSTERIZATION *****
        if options.recluster_outliers:
            with clui.fbm("Performing reclusterization of outliers", cont=False):
                '''
                potentially_recursive_clusterization(config=config,
                                                 clusterization_name=config.recluster_name(),
                                                 inlets_object=inls,
                                                 spaths=spaths,
                                                 traj_reader=traj_reader,
                                                 message='reclusterization',
                                                 max_level=max_level)
                '''
                clustering_function = get_clustering_method(rcoptions, config)
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

        # TODO: Move it after master paths!
        # ***** ADD PASSING PATHS TO CLUSTERS *****
        if options.exclude_passing_in_clusterization and options.add_passing_to_clusters:
            with clui.fbm("Adding passing paths inlets to clusters",cont=False):
                # passing paths were excluded and they are meant to be added
                # one need loop over clusters and then all passing paths have to checked
                # it is assumed taht adding method is barber
                abo = config.get_stage_options(2)  # ab from stage III
                ab_options = get_auto_barber_options(abo)
                ab_options.update(get_auto_barber_options(options))
                # get single paths only (no passing paths)
                spaths_single = [sp for sp in spaths if sp.is_single()]
                # ids of passing paths
                spaths_passing_ids = [nr for nr in xrange(len(spaths)) if spaths[nr].is_passing()]
                if len(spaths_passing_ids) == 0:
                    clui.message("No passing paths to add.")
                else:
                    # loop over passing paths, add to inlets
                    inls.passing = True
                    passing_inlets_ids = []
                    for passing_id in spaths_passing_ids:
                        passing_inlets_ids.extend(inls.extend_inlets(spaths[passing_id]))
                    # loop over clusters
                    for cluster in inls.clusters_list:
                        added_to_cluster = 0
                        if cluster == 0: continue
                        clui.message("Current cluster: %d." % cluster)
                        # sps = inls.lim2clusters(cluster).limspaths2(spaths_single)
                        # chull = inls.lim2clusters(cluster).get_chull()
                        wtc = WhereToCut(inlets=inls.lim2clusters(cluster), **ab_options)
                        wtc.cut_thyself()
                        pbar = clui.SimpleProgressBar(len(passing_inlets_ids),"Loop over available passing paths inlets:")
                        for passing_inlet_nr in range(len(passing_inlets_ids))[::-1]:
                            inlet = inls.inlets_list[passing_inlets_ids[passing_inlet_nr]]
                            sphere = wtc.inlet2sphere(inlet)
                            if sphere is not None:
                                # if True:
                                if wtc.is_overlaping_with_cloud(sphere):
                                    # if chull.point_within(inlet.coords):
                                    # add this inlet to cluster!
                                    inls.clusters[passing_inlets_ids[passing_inlet_nr]] = cluster
                                    added_to_cluster += 1
                                    passing_inlets_ids.pop(passing_inlet_nr)
                            pbar.next()
                        if added_to_cluster:
                            inls.add_message_wrapper(message='+%d passing' % added_to_cluster, toleaf=cluster)
                        pbar.finish()
                    if len(passing_inlets_ids):
                        inls.add_message_wrapper(message='+%d passing' % len(passing_inlets_ids), toleaf=0)

        clui.message('Clustering history:')
        clui.message(clui.print_simple_tree(inls.tree, prefix='').rstrip())

        with clui.fbm("Calculating cluster types"):
            ctypes = inls.spaths2ctypes(spaths)

        # now, there is something to do with ctypes!
        # we can create master paths!
        # but only if user wants this
        master_paths = {}
        master_paths_smooth = {}
        if options.create_master_paths:
            with clui.fbm("Master paths calculations", cont=False):
                smooth = get_smooth_method(soptions)  # this have to preceed GCS
                if GCS.cachedir or GCS.cachemem:
                    pbar = clui.pbar(len(spaths)*2, mess='Building coords cache')
                    [sp.get_coords(smooth=None) for sp in spaths if pbar.next() is None and not isinstance(sp, PassingPath)]
                    [sp.get_coords(smooth=smooth) for sp in spaths if pbar.next() is None and not isinstance(sp, PassingPath)]
                    pbar.finish()
                    use_threads = optimal_threads.threads_count
                else:
                    logger.warning(
                        "Master paths calculation without cache-dir or cache-mem option can be EXTREMELY slow.")
                    use_threads = 1
                with clui.fbm("Creating master paths for cluster types", cont=False):
                    ctypes_generic = [ct.generic for ct in ctypes]
                    ctypes_generic_list = sorted(list(set(ctypes_generic)))

                    pbar = clui.pbar(len([None for sp in spaths if not isinstance(sp, PassingPath)]) * 2)
                    for nr, ct in enumerate(ctypes_generic_list):
                        logger.debug('CType %s (%d)' % (str(ct), nr))
                        sps = lind(spaths, what2what(ctypes_generic, [ct]))
                        # no passing paths are allowed
                        sps = [sp for sp in sps if not isinstance(sp, PassingPath)] # no PassingPaths!
                        if not len(sps):
                            logger.debug(
                                'CType %s (%d), no single paths found, MasterPath calculation skipped.' % (
                                str(ct), nr,))
                            continue
                        logger.debug('CType %s (%d), number of spaths %d' % (str(ct), nr, len(sps)))
                        # print len(sps),ct
                        ctspc = CTypeSpathsCollection(spaths=sps, ctype=ct, pbar=pbar,
                                                      threads=use_threads)
                        master_paths.update({ct: ctspc.get_master_path(resid=(0, nr))})
                        master_paths_smooth.update({ct: ctspc.get_master_path(resid=(0, nr), smooth=smooth)})
                        del ctspc
                    pbar.finish()

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
    if len(template) != len(line):
        pass
    # detect nan and int problem
    for nr, (t, l) in enumerate(zip(template, line)):
        if is_number(l):
            if np.isnan(l):
                if 'd' in t:
                    template[nr] = t.replace('d', 's')
                    line[nr] = 'nan'
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


def spath_name_header():
    return ['RES'], ['%4s']


def add_path_id_head(gen):
    sph, splt = zip(spath_id_header(), spath_name_header())
    sph = [e[0] for e in sph]
    splt = [e[0] for e in splt]

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
            line = [spath.id, spath.id.name] + line
        return line

    return patched


def size_header():
    return ['Size', 'Size%'], ['%7d', '%6.2f']


def add_size_head(gen):
    sph, splt = size_header()

    @wraps(gen)
    def patched(*args, **kwargs):
        h, lt = gen(*args, **kwargs)
        return sph + h, splt + lt

    return patched


def add_size(gen):
    @wraps(gen)
    def patched(spaths, add_size=True, add_size_p100=None, *args, **kwargs):
        line = gen(spaths, *args, **kwargs)
        if add_size_p100 is not None:
            line = [len(spaths) / float(add_size_p100) * 100] + line
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

    # TODO: Change it in such a way that it cooperates well with debug-file option.
    def __init__(self, fileoption, line_nr=False):
        self.output2stderr = False
        if fileoption:
            self.filehandle = open(fileoption, 'w')
            # self.output2stderr = True
        else:
            self.filehandle = sys.stdout
        self.line_nr = line_nr

    def __call__(self, info2print, nr=None):
        if self.line_nr and nr is not None:
            info2print = (self.nr_template % nr) + info2print
        if self.output2stderr:
            clui.message(info2print)
        print >> self.filehandle, info2print

    def sep(self):
        self(asep())

    def thead(self, info2print):
        self(clui.thead(info2print))

    def tend(self, info2print):
        self(clui.tsep(info2print))

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
    if not isinstance(spath, PassingPath):
        line.extend(map(len, (spath.path_in, spath.path_object, spath.path_out)))
    else:
        line += [float('nan')] * 3
    line.append(spath.ends)
    return line


################

@add_path_id_head
def spath_lenght_total_info_header(total=None):
    header = 'InpL ObjL OutL'.split()
    if total:
        header = ['TotL'] + header
    line_template = ['%9.1f'] * len(header)
    return header, line_template


@add_path_id
def spath_lenght_total_info(spath, totalonly=False, total=False):
    # total and totalonly are internal flags
    # to calculate: total len, in len, obj len, out len call:
    # total=True (this will calculate total)
    # to skip calculation of total call:
    # total=False, totalonly=False (this will calculate in, obj, and out lens)
    # to calculate total len call:
    # total=False, totalonly=True
    line = []
    if not total:
        if not totalonly:
            for t in traces.midpoints(spath.coords):
                if len(t) > 1:  # traces.length_step_std requires at least 2 points
                    line.append(traces.length_step_std(t)[0])
                else:
                    line.append(float('nan'))
            return line
        else:
            for t in traces.midpoints((spath.coords_cont,)):
                if len(t) > 1:  # traces.length_step_std requires at least 2 points
                    line.append(traces.length_step_std(t)[0])
                else:
                    line.append(float('nan'))
            return line
    line += spath_lenght_total_info(spath, add_id=False, total=False, totalonly=True)
    if not isinstance(spath, PassingPath):
        line += spath_lenght_total_info(spath, add_id=False, total=False, totalonly=False)
    else:
        line += [float('nan')] * 3
    return line


@add_path_id
def spath_frames_total_info(spath, totalonly=False, total=False):
    # total and totalonly are internal flags
    # to calculate: total len, in len, obj len, out len call:
    # total=True (this will calculate total)
    # to skip calculation of total call:
    # total=False, totalonly=False (this will calculate in, obj, and out lens)
    # to calculate total len call:
    # total=False, totalonly=True
    line = []
    if not total:
        if not totalonly:
            for t in spath.coords:
                line.append(len(t))
            return line
        else:
            for t in (spath.coords_cont,):
                line.append(len(t))
            return line
    line += spath_frames_total_info(spath, add_id=False, total=False, totalonly=True)
    if not isinstance(spath, PassingPath):
        line += spath_frames_total_info(spath, add_id=False, total=False, totalonly=False)
    else:
        line += [float('nan')] * 3
    return line


################

@add_path_id_head
def spath_steps_info_header(total=None):
    header = 'InpS InpStdS ObjS ObjStdS OutS OutStdS'.split()
    if total:
        header = ['TotS', 'TotStdS'] + header
    line_template = ['%8.2f', '%8.3f'] * (len(header) / 2)
    return header, line_template


@add_path_id
def spath_steps_info(spath, total=None):
    line = []
    if not total:
        for t in traces.midpoints(spath.coords):
            if len(t) > 0:
                line.extend(traces.length_step_std(t)[1:])
            else:
                line.extend([float('nan'), float('nan')])
        return line
    line += spath_steps_info(spath, add_id=False, total=False)
    if not isinstance(spath, PassingPath):
        t = traces.midpoints((spath.coords_cont,)).next()
        if len(t) > 0:
            line = list(traces.length_step_std(t)[1:]) + line
        else:
            line = [float('nan'), float('nan')] + line
    else:
        line += [float('nan'), float('nan')] * 3
    return line


@add_path_id
def spath_frames_info(spath, total=None):
    line = []
    if not total:
        for t in spath.coords:
            line.append(len(t))
        return line
    line += spath_steps_info(spath, add_id=False, total=False)
    if not isinstance(spath, PassingPath):
        t = spath.coords_cont
        line = [len(t)] + line
    else:
        line = [float('nan')] + line
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
def spath_full_info_header(total=None):
    header = []
    line_template = []
    for h, lt in (spath_basic_info_header(add_id=False),
                  spath_lenght_total_info_header(add_id=False, total=total),
                  spath_steps_info_header(add_id=False, total=total),
                  spath_ctype_header(add_id=False)):
        header += h
        line_template += lt
    return header, line_template


@add_path_id
def spath_full_info(spath, ctype=None, total=None):
    line = []
    for l in (spath_basic_info(spath, add_id=False),
              spath_lenght_total_info(spath, add_id=False, total=total),
              spath_steps_info(spath, add_id=False, total=total),
              spath_ctype(spath, ctype=ctype, add_id=False)):
        line += l
    return line


################################################################################

@add_size_head
def spaths_lenght_total_header():
    header = 'Tot TotStd Inp InpStd Obj ObjStd Out OutStd'.split()
    line_template = ['%9.1f', '%9.2f'] * (len(header) / 2)
    return header, line_template


@add_size
def spaths_length_total(spaths):
    line = []
    d4s = []
    for sp in spaths:
        if not isinstance(sp, PassingPath):
            d4s.append(spath_lenght_total_info(sp, add_id=False))
    d4s = np.array(d4s)
    if d4s.size:
        line.extend(np.mean(d4s, 0))
        line.extend(np.std(d4s, 0))
    else:
        line = [float('nan')] * 6
    # total
    d4s = []
    for sp in spaths:
        d4s.append(spath_lenght_total_info(sp, totalonly=True, add_id=False))
    d4s = np.array(d4s)
    if d4s.size:
        line.extend(np.mean(d4s, 0))
        line.extend(np.std(d4s, 0))
    else:
        line += [float('nan')] * 4

    return [line[6], line[7], line[0], line[3], line[1], line[4], line[2], line[5]]


@add_size
def spaths_frames_total(spaths):
    line = []
    d4s = []
    for sp in spaths:
        if not isinstance(sp, PassingPath):
            d4s.append(spath_frames_total_info(sp, add_id=False))
    d4s = np.array(d4s)
    if d4s.size:
        line.extend(np.mean(d4s, 0))
        line.extend(np.std(d4s, 0))
    else:
        line = [float('nan')] * 6
    # total
    d4s = []
    for sp in spaths:
        d4s.append(spath_frames_total_info(sp, totalonly=True, add_id=False))
    d4s = np.array(d4s)
    if d4s.size:
        line.extend(np.mean(d4s, 0))
        line.extend(np.std(d4s, 0))
    else:
        line += [float('nan')] * 4

    return [line[6], line[7], line[0], line[3], line[1], line[4], line[2], line[5]]


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


################################################################################

@add_ctype_id_head
def ctypes_spaths_info_header():
    header, line_template = spaths_lenght_total_header()
    return header, line_template


@add_ctype_id
def ctypes_spaths_info(ctype, spaths, show='len', add_size_p100=None):
    # show could be len or frames
    line = []
    if show == 'len':
        line += spaths_length_total(spaths, add_size_p100=add_size_p100)
    if show == 'frames':
        line += spaths_frames_total(spaths, add_size_p100=add_size_p100)
    return line


################################################################################


@add_cluster_id_head
def clusters_stats_prob_header():
    header = 'IN-OUT diff N IN-OUT_prob diff_prob N_prob'.split()
    line_template = ['%8d'] * (len(header) / 2) + ['%12.2f'] * (len(header) / 2)
    # header += 'IN_len OUT_len Both_len'.split()
    # line_template += ['%9.1f'] * 3
    return header, line_template


@add_cluster_id
def clusters_stats_prob(cluster, sp_ct):
    # calculates probabilities of some events for cluster
    # X:X transition - io
    # X:? and ?:X transition - d
    # X:N and N:X transition - N
    line = []
    io, d, N = 0, 0, 0
    in_len, out_len, tot_len = 0., 0., 0.
    in_n, out_n, tot_n = 0, 0, 0
    for sp, ct in sp_ct:
        ct = ct.generic.clusters
        assert cluster in ct
        if cluster == ct[0] and cluster == ct[1]:
            io += 1
        elif None in ct:
            N += 1
        else:
            d += 1
    line += [io, d, N]
    summa = float(sum([io, d, N]))
    line += map(lambda x: x / summa if summa else float('nan'), [io, d, N])
    return line


@add_cluster_id_head
def clusters_stats_len_header():
    header = 'X->Obj Obj->X p-value X->ObjMin X->ObjMinID Obj->XMin Obj->XMinID'.split()
    line_template = (['%9.1f'] * 2) + ['%9.4f'] + ['%9.1f', '%11s'] * 2
    return header, line_template


@add_cluster_id
def clusters_stats_len(cluster, sp_ct):
    line = []
    in_len, out_len = [], []
    in_len_min, out_len_min = float('inf'), float('inf')
    in_len_min_id, out_len_min_id = None, None
    for sp, ct in sp_ct:
        ct = ct.generic.clusters
        assert cluster in ct
        lens = spath_lenght_total_info(sp, add_id=False, total=False, totalonly=False)
        # tot,in,obj,out
        if cluster == ct[0]:
            if not np.isnan(lens[0]):
                in_len.append(lens[0])
                if lens[0] < in_len_min:
                    in_len_min = lens[0]
                    in_len_min_id = str(sp.id)
        if cluster == ct[1]:
            if not np.isnan(lens[-1]):
                out_len.append(lens[-1])
                if lens[-1] < out_len_min:
                    out_len_min = lens[-1]
                    out_len_min_id = str(sp.id)

    if len(in_len):
        line.append(np.mean(in_len))
    else:
        line.append(float('nan'))
    if len(out_len):
        line.append(np.mean(out_len))
    else:
        line.append(float('nan'))
    if len(in_len) > 1 and len(out_len) > 1:
        line.append(ttest_ind(in_len, out_len)[-1])  # this is supposed to return p-value
    else:
        line.append(float('nan'))

    line += [in_len_min, in_len_min_id, out_len_min, out_len_min_id]

    return line


@add_cluster_id_head
def clusters_stats_steps_header():
    header = 'X->Obj Obj->X p-value X->ObjMin X->ObjMinID Obj->XMin Obj->XMinID'.split()
    line_template = (['%9.1f'] * 2) + ['%9.4f'] + ['%9.1f', '%11s'] * 2
    return header, line_template


@add_cluster_id
def clusters_stats_steps(cluster, sp_ct):
    line = []
    in_len, out_len = [], []
    in_len_min, out_len_min = float('inf'), float('inf')
    in_len_min_id, out_len_min_id = None, None
    for sp, ct in sp_ct:
        ct = ct.generic.clusters
        assert cluster in ct
        lens = spath_frames_info(sp, add_id=False, total=False)
        # tot,in,obj,out
        if cluster == ct[0]:
            if not np.isnan(lens[0]):
                in_len.append(lens[0])
                if lens[0] < in_len_min:
                    in_len_min = lens[0]
                    in_len_min_id = str(sp.id)
        if cluster == ct[1]:
            if not np.isnan(lens[-1]):
                out_len.append(lens[-1])
                if lens[-1] < out_len_min:
                    out_len_min = lens[-1]
                    out_len_min_id = str(sp.id)

    if len(in_len):
        line.append(np.mean(in_len))
    else:
        line.append(float('nan'))
    if len(out_len):
        line.append(np.mean(out_len))
    else:
        line.append(float('nan'))
    if len(in_len) > 1 and len(out_len) > 1:
        line.append(ttest_ind(in_len, out_len)[-1])  # this is supposed to return p-value
    else:
        line.append(float('nan'))

    line += [in_len_min, in_len_min_id, out_len_min, out_len_min_id]

    return line


################################################################################


# analysis
def stage_V_run(config, options,
                spaths=None,
                paths=None,
                inls=None,
                ctypes=None,
                reader=None,
                **kwargs):
    # file handle?
    head_nr = False
    line_nr = head_nr
    pa = PrintAnalysis(options.save, line_nr=line_nr)

    if options.save:
        clui.message('Using user provided file (%s), and' % options.save)
        clui.message('for histograms data file (%s).' % (options.save + '.csv'))
    else:
        clui.message('Using standard output.')
        clui.message(sep())
        clui.message('')

    ############
    pa.sep()
    pa('Aqua-Duct analysis')
    pa(clui.get_str_timestamp())

    ############
    if options.dump_config:
        pa.sep()
        pa.under('Configuration file name: %s' % config.config_filename)
        pa(os.linesep.join(config.dump_config()))

    ############
    pa.sep()
    pa("Frames window: %d:%d step %d" % (Reader.window.start,
                                         Reader.window.stop,
                                         Reader.window.step))

    ############
    traced_names = tuple(sorted(list(set([sp.id.name for sp in spaths]))))
    # traced_names = ['all'] + traced_names

    pa.sep()
    pa("Names of traced molecules: %s" % (' '.join(traced_names)))

    # spaths_types = ['all']
    spaths_types = []
    for sp in spaths:
        if isinstance(sp, PassingPath):
            if not PassingPath in spaths_types:
                spaths_types.append(PassingPath)
        else:
            if not SinglePath in spaths_types:
                spaths_types.append(SinglePath)
    spaths_types = tuple(spaths_types)

    ############

    def iter_over_tn():
        yield traced_names, ''
        if len(traced_names) > 1:
            for _tname in traced_names:
                yield (_tname,), " of %s" % _tname

    def iter_over_tnspt():
        for _tname, _message in iter_over_tn():
            yield _tname, spaths_types, _message
            if len(spaths_types) > 1:
                for _sptype in spaths_types:
                    if _sptype == PassingPath:
                        _sptype_name = 'passing paths'
                    elif _sptype == SinglePath:
                        _sptype_name = 'object  paths'
                    yield _tname, (_sptype,), _message + (" of %s" % _sptype_name)

    ############

    pa.sep()
    for tname, message in iter_over_tn():
        pa("Number of traceable residues%s: %d" %
           (message,
            len([None for p in paths if p.name in tname])))

    for tname, sptype, message in iter_over_tnspt():
        pa("Number of separate paths%s: %d" %
           (message,
            len([None for sp in spaths if (sp.id.name in tname) and
                 (isinstance(sp, sptype))])))

    ############
    pa.sep()
    for tname, sptype, message in iter_over_tnspt():
        pa("Number of inlets%s: %d" % (
            message, inls.lim2spaths([sp for sp in spaths if isinstance(sp, sptype)]).lim2rnames(tname).size))

    no_of_clusters = len(inls.clusters_list) - {True: 1, False: 0}[0 in inls.clusters_list]  # minus outliers, if any
    pa("Number of clusters: %d" % no_of_clusters)
    pa("Outliers: %s" % ({True: 'yes', False: 'no'}[0 in inls.clusters_list]))

    ############
    pa.sep()
    pa('Clustering history:')
    pa(clui.print_simple_tree(inls.tree, prefix='').rstrip())

    ############
    pa.sep()
    header_line, line_template = get_header_line_and_line_template(clusters_inlets_header(), head_nr=head_nr)
    for tname, sptype, message in iter_over_tnspt():
        pa("Clusters summary - inlets%s" % message)
        pa.thead(header_line)
        for nr, cl in enumerate(inls.clusters_list):
            inls_lim = inls.lim2spaths([sp for sp in spaths if isinstance(sp, sptype)]).lim2rnames(tname).lim2clusters(
                cl)
            pa(make_line(line_template, clusters_inlets(cl, inls_lim)), nr=nr)
        pa.tend(header_line)

    ############
    ############
    ############
    pa.sep()

    for tname, sptype, message in iter_over_tnspt():
        header_line, line_template = get_header_line_and_line_template(clusters_stats_prob_header(), head_nr=head_nr)
        pa("Clusters statistics (of paths%s) probabilities of transfers" % message)
        pa.thead(header_line)
        for nr, cl in enumerate(inls.clusters_list):
            sp_ct_lim = ((sp, ct) for sp, ct in zip(spaths, ctypes) if
                         cl in ct.clusters and isinstance(sp, sptype) and sp.id.name in tname)
            pa(make_line(line_template, clusters_stats_prob(cl, sp_ct_lim)), nr=nr)
        pa.tend(header_line)

        header_line, line_template = get_header_line_and_line_template(clusters_stats_len_header(), head_nr=head_nr)
        pa("Clusters statistics (of paths%s) mean lengths of transfers" % message)
        pa.thead(header_line)
        for nr, cl in enumerate(inls.clusters_list):
            sp_ct_lim = ((sp, ct) for sp, ct in zip(spaths, ctypes) if
                         cl in ct.clusters and isinstance(sp, sptype) and sp.id.name in tname)
            pa(make_line(line_template, clusters_stats_len(cl, sp_ct_lim)), nr=nr)
        pa.tend(header_line)

        header_line, line_template = get_header_line_and_line_template(clusters_stats_steps_header(), head_nr=head_nr)
        pa("Clusters statistics (of paths%s) mean frames numbers of transfers" % message)
        pa.thead(header_line)
        for nr, cl in enumerate(inls.clusters_list):
            sp_ct_lim = ((sp, ct) for sp, ct in zip(spaths, ctypes) if
                         cl in ct.clusters and isinstance(sp, sptype) and sp.id.name in tname)
            pa(make_line(line_template, clusters_stats_steps(cl, sp_ct_lim)), nr=nr)
        pa.tend(header_line)

    ############
    pa.sep()
    header_line, line_template = get_header_line_and_line_template(ctypes_spaths_info_header(), head_nr=head_nr)

    ctypes_generic = [ct.generic for ct in ctypes]
    ctypes_generic_list = sorted(list(set(ctypes_generic)))

    # sorted by ctype
    ctypes_size = []
    for nr, ct in enumerate(ctypes_generic_list):
        sps = lind(spaths, what2what(ctypes_generic, [ct]))
        ctypes_size.append(len(sps))
    ctypes_generic_list = [ctypes_generic_list[i] for i in np.argsort(ctypes_size)[::-1]]

    for tname, sptype, message in iter_over_tnspt():
        pa("Separate paths clusters types summary - mean lengths of paths%s" % message)
        pa.thead(header_line)

        total_size = len([sp for sp in spaths if sp.id.name in tname and isinstance(sp, sptype)])
        for nr, ct in enumerate(ctypes_generic_list):
            sps = lind(spaths, what2what(ctypes_generic, [ct]))
            sps = [sp for sp in sps if sp.id.name in tname and isinstance(sp, sptype)]
            # ctypes_size.append(len(sps))
            if len(sps) > 0:
                pa(make_line(line_template, ctypes_spaths_info(ct, sps, add_size_p100=total_size, show="len")), nr=nr)
        pa.tend(header_line)
        pa("Separate paths clusters types summary - mean number of frames of paths%s" % message)
        pa.thead(header_line)
        for nr, ct in enumerate(ctypes_generic_list):
            sps = lind(spaths, what2what(ctypes_generic, [ct]))
            sps = [sp for sp in sps if sp.id.name in tname and isinstance(sp, sptype)]
            # ctypes_size.append(len(sps))
            if len(sps) > 0:
                pa(make_line(line_template, ctypes_spaths_info(ct, sps, add_size_p100=total_size, show="frames")),
                   nr=nr)
        pa.tend(header_line)

    ############
    pa.sep()
    pa("List of separate paths and properties")
    header_line, line_template = get_header_line_and_line_template(spath_full_info_header(total=True), head_nr=head_nr)
    pa.thead(header_line)
    for nr, (sp, ctype) in enumerate(izip_longest(spaths, ctypes, fillvalue=None)):
        if ctype is not None:
            ctype = ctype.generic
        pa(make_line(line_template, spath_full_info(sp, ctype=ctype, total=True)), nr=nr)
    pa.tend(header_line)

    ############
    # additional analysis

    def iter_over_tn():
        yield traced_names, 'amol'
        if len(traced_names) > 1:
            for _tname in traced_names:
                yield (_tname,), "%s" % _tname

    def iter_over_part():
        for _part in 'walk in object out'.split() + ['in out'.split()]:
            if isinstance(_part, list):
                _message = '_'.join(_part)
            else:
                _message = _part
                _part = [_part]
            yield _part, _message

    def iter_over_spt():
        yield spaths_types, 'apaths'
        if len(spaths_types) > 1:
            for _sptype in spaths_types:
                if _sptype == PassingPath:
                    _sptype_name = 'passing'
                elif _sptype == SinglePath:
                    _sptype_name = 'object'
                yield (_sptype,), _sptype_name

    def iter_over_c():
        yield inls.clusters_list + [None], 'aclusts'
        if len(inls.clusters_list + [None]) > 1:
            for _cluster in inls.clusters_list + [None]:
                yield [_cluster], str(_cluster)

    def iter_over_ct():
        yield ctypes_generic_list, 'actypes'
        if len(ctypes_generic_list) > 1:
            for _cluster in ctypes_generic_list:
                yield (_cluster,), str(_cluster)

    def iter_over_all():
        for _tname, _m_tname in iter_over_tn():
            for _sptype, _m_sptype in iter_over_spt():
                for _cluster, _m_cluster in iter_over_c():
                    for _part, _m_part in iter_over_part():
                        _message = '_'.join((_m_tname,
                                             _m_sptype,
                                             _m_cluster,
                                             _m_part))
                        yield _tname, _sptype, _cluster, _part, _message
                for _ctype, _m_ctype in iter_over_ct():
                    for _part, _m_part in iter_over_part():
                        _message = '_'.join((_m_tname,
                                             _m_sptype,
                                             _m_ctype,
                                             _m_part))
                        yield _tname, _sptype, _ctype, _part, _message

    # histograms
    # calculate old max_frame
    max_frame = Reader.number_of_frames(onelayer=True)
    header = [column[-1] for column in iter_over_all()]
    fmt = ['%u'] * len(header)
    h = np.zeros((max_frame, len(header)))
    # loop over spaths
    pbar = clui.pbar(maxval=len(spaths),
                     mess='Calculating histograms')
    # loop over paths and ctypes
    for sp, ct in zip(spaths, ctypes):
        # loop over columns
        # traced names, paths types, clusters cluster types, part of paths, column name
        for tname, sptype, c_ct, part, col_name in iter_over_all():
            # check if column fits to the requirements
            if not sp.id.name in tname: continue
            if not isinstance(sp, sptype): continue
            # c or ct
            it_is_ct = False
            if isinstance(c_ct[0], InletClusterGenericType):
                it_is_ct = True
                if not ct.generic in c_ct: continue
            else:
                if len(union_full(ct.generic.clusters, c_ct)) == 0: continue
            # this is more than that...
            # if c_ct is not InletClusterGenericType then:
            # if 'in' in part only incoming paths in ct are used
            # if 'out' in part only outgoing paths in ct are used
            col_index = header.index(col_name)
            if 'walk' in part:
                h[sp.paths_cont, col_index] += 1
            if isinstance(sp, PassingPath): continue
            if 'in' in part:
                if not it_is_ct:
                    if not ct.input in c_ct:
                        continue
                h[sp.path_in, col_index] += 1
            if 'object' in part:
                h[sp.path_object, col_index] += 1
            if 'out' in part:
                if not it_is_ct:
                    if not ct.output in c_ct:
                        continue
                h[sp.path_out, col_index] += 1
        pbar.next()
    pbar.finish()

    # scope/object size?
    if options.calculate_scope_object_size:
        scope_size = []
        object_size = []
        # now, the problem is in the scope and object definition.
        pbar = clui.pbar(maxval=Reader.number_of_frames(), mess='Calculating scope and object sizes')

        for number, traj_reader in Reader.iterate(number=True):
            scope_size.append([])
            object_size.append([])
            for frame in traj_reader.iterate_over_frames():
                scope = traj_reader.parse_selection(options.scope_chull)
                ch = scope.chull()
                scope_size[-1].append((ch.area, ch.volume))
                res = traj_reader.parse_selection(options.object_chull)
                ch = res.chull()
                object_size[-1].append((ch.area, ch.volume))
                pbar.next()
            header += map(lambda s: '%s_%d' % (s, number),
                          ['scope_area', 'scope_volume', 'object_area', 'object_volume'])
            fmt += ['%0.3f', '%0.2f'] * 2
        for s_s, o_s in zip(scope_size, object_size):
            h = np.hstack((h, s_s, o_s))
        pbar.finish()
    # add frame column?
    frame_col = np.array([range(max_frame)]).T
    h = np.hstack((frame_col, h))
    header = ['frame'] + header
    fmt = ['%u'] + fmt
    # save???
    if options.save:
        h_fname = options.save + '.csv'
    else:
        import cStringIO as StringIO
        h_fname = StringIO.StringIO()
    np.savetxt(h_fname, h,
               fmt=fmt,
               delimiter=',',
               header=','.join(header))
    if not options.save:
        print h_fname.getvalue()

    return {'hist': h, 'header': header}


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
            # FIXME: if path(s) is empty this will probably produce CGO_END without CGO_BEGIN
            spp.paths_trace(sp, name=name + '_in' + name_separate, plot_walk=False, plot_object=False, plot_out=False,
                            state=state, smooth=smooth)
            spp.paths_trace(sp, name=name + '_obj' + name_separate, plot_walk=False, plot_in=False, plot_out=False,
                            state=state, smooth=smooth)
            spp.paths_trace(sp, name=name + '_out' + name_separate, plot_walk=False, plot_in=False, plot_object=False,
                            state=state, smooth=smooth)
            spp.paths_trace(sp, name=name + '_walk' + name_separate, plot_in=False, plot_object=False, plot_out=False,
                            state=state, smooth=smooth)
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
                 spaths=None,
                 inls=None,
                 ctypes=None,
                 master_paths=None,
                 master_paths_smooth=None,
                 **kwargs):
    from aquaduct.visual.pymol_connector import ConnectToPymol, SinglePathPlotter
    # from aquaduct.visual.pymol_connector import cmd as pymol_cmd
    from aquaduct.visual.helpers import ColorMapDistMap

    soptions = config.get_smooth_options()
    smooth = get_smooth_method(soptions)

    # start pymol
    with clui.fbm("Starting PyMOL connection", cont=False):
        pymol_connector = ConnectToPymol()
        if is_pymol_connector_script(options.save):
            pymol_connector.init_script(options.save)
        else:
            pymol_connector.init_pymol()

        if options.simply_smooths:
            spp = SinglePathPlotter(pymol_connector, linearize=get_linearize_method(options.simply_smooths))
        else:
            spp = SinglePathPlotter(pymol_connector, linearize=None)

    ctypes_generic = [ct.generic for ct in ctypes]
    ctypes_generic_list = sorted(list(set(ctypes_generic)))

    if options.show_molecule:
        molecule_name = ''
        with clui.fbm("Molecule"):
            for nr, traj_reader in enumerate(Reader.iterate()):
                # mda_ppr = mda.core.flags["permissive_pdb_reader"]
                # mda.core.flags["permissive_pdb_reader"] = False #mda16 it is porbably always True
                frames_to_show = range2int(options.show_molecule_frames)
                pdbfile = traj_reader.dump_frames(frames_to_show, selection=options.show_molecule)
                pymol_connector.load_pdb('molecule%d' % nr, pdbfile)
                if len(molecule_name) == 0:
                    molecule_name = 'molecule%d' % nr
                os.unlink(pdbfile)
                # it would be nice to plot convexhull
    if options.show_scope_chull:
        with clui.fbm("Convexhull"):
            for nr, traj_reader in enumerate(Reader.iterate()):
                frames_to_show = range2int(options.show_scope_chull_frames)
                for frame in frames_to_show:
                    traj_reader.set_frame(frame)
                    scope = traj_reader.parse_selection(options.show_scope_chull)
                    chull = scope.chull()
                    spp.convexhull(chull, name='scope_shape%d' % nr, state=frame + 1)

    if options.show_object_chull:
        with clui.fbm("Object shape"):
            for nr, traj_reader in enumerate(Reader.iterate()):
                frames_to_show = range2int(options.show_scope_chull_frames)
                for frame in frames_to_show:
                    traj_reader.set_frame(frame)
                    object_shape = traj_reader.parse_selection(options.show_object_chull)
                    chull = object_shape.chull()
                    spp.convexhull(chull, name='object_shape%d' % nr, color=np.array([255, 153, 0]) / 255.,
                                   state=frame + 1)  # orange

    if options.inlets_clusters:
        with clui.fbm("Clusters"):
            # TODO: require stage V for that?
            no_of_clusters = len(inls.clusters_list)  # total, including outliers
            cmap = ColorMapDistMap()
            for c in inls.clusters_list:
                # coords for current cluster
                ics = inls.lim2clusters(c).coords
                if c == 0:
                    c_name = 'out'
                else:
                    c_name = str(int(c))
                spp.scatter(ics, color=cmap(c), name="cluster_%s" % c_name)
                if False:  # TODO: This does not work any more in that way. Rewrite it or remove it
                    radii = inls.lim2clusters(c).radii
                    if len(radii) > 0:
                        spp.scatter(ics, color=cmap(c), radius=radii, name="cluster_radii_%s" % c_name)

    if options.ctypes_raw:
        with clui.fbm("CTypes raw"):
            for nr, ct in enumerate(ctypes_generic_list):
                clui.message(str(ct), cont=True)
                sps = lind(spaths, what2what(ctypes_generic, [ct]))
                plot_spaths_traces(sps, name=str(ct) + '_raw', split=False, spp=spp)
                if ct in master_paths:
                    if master_paths[ct] is not None:
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
        pymol_connector.orient_on(molecule_name)

    if is_pymol_connector_session(options.save):
        with clui.fbm("Saving session (%s)" % options.save):
            clui.message("")  # new line
            pbar = clui.pbar(len(spaths))
            # FIXME: Loop over states is not required if there is no object with many states.
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


__all__ = '''clui
optimal_threads
aquaduct_version_nice
ValveConfig
valve_begin
valve_load_config
valve_exec_stage
stage_I_run
stage_II_run
stage_III_run
stage_IV_run
stage_V_run
stage_VI_run
valve_end
'''.split()
