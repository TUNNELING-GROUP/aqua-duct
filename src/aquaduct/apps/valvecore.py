# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2016-2017  Tomasz Magdziarz, Alicja Płuciennik, Michał Stolarczyk <info@aquaduct.pl>
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


from aquaduct import logger, logger_name

################################################################################

import multiprocessing as mp
import ConfigParser
import cPickle as pickle
import StringIO
import copy
import gzip
import os
import re
import operator
import shlex
import sys

import numpy as np
import MDAnalysis as mda
import roman  # TODO: remove this dependency!

from collections import namedtuple, OrderedDict  # TODO: check if OrderedDict is REALLY used
from functools import wraps
from itertools import izip_longest
from keyword import iskeyword
from distutils.version import LooseVersion

from scipy.spatial.distance import cdist
# If scipy is relatively old and numpy is relatively new this triggers warning on oldnumeric module deprecation.
# This warning emerges if MDAnalysis is imported and then scipy. Observed under Ubuntu 14.04.

from aquaduct import greetings as greetings_aquaduct
from aquaduct import version as aquaduct_version
from aquaduct import version_nice as aquaduct_version_nice
from aquaduct.geom import traces
from aquaduct.geom.cluster import PerformClustering, DBSCAN, AffinityPropagation, MeanShift, KMeans, Birch, BarberCluster, get_required_params
from aquaduct.geom.cluster import  AVAILABLE_METHODS as available_clusterization_methods

from aquaduct.geom.convexhull import is_point_within_convexhull
from aquaduct.geom.master import CTypeSpathsCollection
from aquaduct.geom.smooth import WindowSmooth, MaxStepSmooth, WindowOverMaxStepSmooth, ActiveWindowSmooth, \
    ActiveWindowOverMaxStepSmooth, DistanceWindowSmooth, DistanceWindowOverMaxStepSmooth, SavgolSmooth
from aquaduct.traj.dumps import TmpDumpWriterOfMDA
from aquaduct.traj.inlets import Inlets, InletTypeCodes
from aquaduct.traj.paths import GenericPaths, yield_single_paths
from aquaduct.traj.reader import ReadViaMDA,atom2vdw_radius
from aquaduct.traj.selections import CompactSelectionMDA, SelectionMDA
from aquaduct.utils import clui
from aquaduct.utils.helpers import range2int, Auto, what2what, lind, is_number
from aquaduct.utils.multip import optimal_threads
from aquaduct.traj.barber import WhereToCut

def version():
    # TODO: remove it
    return 0, 10, 0


def version_nice():
    return '.'.join(map(str, version()))


def version_onenumber():
    return ''.join(map(str, version()))


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
            if opt.replace('_','').isalnum():
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
        #config = RawConfigParserWithComments()

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

        def common_traj_data(section,commented=False):
            for setting in self.common_traj_data_config_names():
                if commented:
                    setting = '#'+setting
                config.set(section, setting, 'None')

        ################
        # global settings
        section = self.global_name()
        config.add_section(section)
        #config.add_comment(section,'global settings')


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
        config.set(section, 'discard_short_paths', '1')

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
        config.set(section, 'show_chull', 'None')
        config.set(section, 'show_chull_frames', '0')
        config.set(section, 'show_object', 'None')
        config.set(section, 'show_object_frames', '0')

        return config

    def load_config(self, filename):
        self.config.read(filename)
        self.config_filename = filename

    def save_config_stream(self, fs):
        self.config.write(fs)

    def save_config(self, filename):
        with open(filename, 'w') as fs:
            self.save_config_stream(fs)

    def get_general_comment(self,section):
        if section == self.global_name():
            return ['Global settings.']
        if section == self.cluster_name():
            out = ['Possible clusterization methods:']
            # get default clusterization method
            defmet = self.config.get(section,'method')
            for method in available_clusterization_methods:
                if method == defmet: continue
                out.append('method = %s' %method)
                params = get_required_params(method)
                if params:
                    sss = ''
                    if len(params) > 1:
                        sss = 's'
                    out.append('Required parameter%s for %s method:' % (sss,method))
                    for param in params:
                        out.append('%s = None' % param)
                    #out.append('')
            return out


    def dump_config(self,dump_template=False):
        skip_list = '''simply_smooths'''.split() # these options are very optional!
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
            output.append('[%s]' % name) # section name
            # general comments
            if dump_template:
                comment = self.get_general_comment(name)
                if comment:
                    output.extend(['# %s' % line for line in comment])
            for key,value in opts._asdict().iteritems(): # loop over options
                if key in skip_list: continue
                # comment scope etc. in stage II if dump_template
                if dump_template and name == self.stage_names(1):
                    if key in self.common_traj_data_config_names():
                        key = '#' + key
                output.append('%s = %s' % (key, value2str(value)))
            if not concise: output.append('')

        # is something missing? another loop over all additional sections
        for miss in self.config.sections():
            if miss in names: continue # skip if it was already dumped in the loop above
            output.append('[%s]' % miss)
            for key in self.config.options(miss): # loop over options in section miss
                if key in skip_list: continue
                value = self.config.get(miss, key)
                output.append('%s = %s' % (key, value2str(value)))
            if not concise: output.append('')
        while not len(output[-1].strip()):
            output.pop()
        return output


################################################################################
# reader helper class

class TrajectoryReader(object):
    def __init__(self, top, trj,frames_window=None):
        assert isinstance(top, (str, unicode)), "Topology file name missing, %s given instead" % str(top)
        assert isinstance(trj, (str, unicode)), "Trajectory file(s) name(s) missing, %s given instead" % str(trj)
        self.top = top
        self.trj = shlex.split(trj)
        self.frames_window = frames_window

    def get(self):
        return ReadViaMDA(self.top, self.trj, window=self.frames_window)

    @property
    def max_frame(self):
        with self.get() as tmp_reader:
            return tmp_reader.number_of_frames - 1  # returns 0-based value


def rebuild_selection(selection, reader):
    return CompactSelectionMDA(selection).toSelectionMDA(reader)


################################################################################
# convex hull helpers
# TODO: Move it to separate module.
# TODO: Following functions are or will be deprecated, remove them as soon as possible.

def CHullCheck(point):
    return CHullCheck.chull.point_within(point)


def CHullCheck_init(args):
    CHullCheck.chull = copy.deepcopy(args[0])


def CHullCheck_pool(chull, threads=optimal_threads.threads_count):
    return mp.Pool(threads, CHullCheck_init, [(chull,)])


def CHullCheck_exec(chull, points, threads=optimal_threads.threads_count):
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
        is_res_in_scope = CHullCheck_exec(chull, res_coords, threads=optimal_threads.threads_count)
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
                         'aquaduct_version': aquaduct_version()}, f)
            # then data to save:
            pickle.dump(data_to_save, f)
            # then other kwargs
            pickle.dump(kwargs, f)


class LoadDumpWrapper(object):
    """This is wrapper for pickled data that provides compatibility
    with earlier versions of Aqua-Duct.

    Conversions in use:

    1) replace 'aquaduct.' by 'aquaduct.'

    """

    def __init__(self, filehandle):
        self.fh = filehandle

    def convert(self, s):
        new_s = s
        new_s = new_s.replace('aqueduct.', 'aquaduct.')
        new_s = new_s.replace('aqueduct_version', 'aquaduct_version')
        return new_s

    def read(self, *args, **kwargs):
        return self.convert(self.fh.read(*args, **kwargs))

    def readline(self, *args, **kwargs):
        return self.convert(self.fh.readline(*args, **kwargs))


def load_dump(filename):
    with clui.fbm('Loading data dump from %s file' % filename):
        with gzip.open(filename, mode='r') as protof:
            f = LoadDumpWrapper(protof)
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


def check_version_compliance(current, loaded, what):
    if current[0] > loaded[0]:
        logger.error('Loaded data has %s major version lower then the application.' % what)
    if current[0] < loaded[0]:
        logger.error('Loaded data has %s major version higher then the application.' % what)
    if current[0] != loaded[0]:
        logger.error('Possible problems with API compliance.')
    if current[1] > loaded[1]:
        logger.warning('Loaded data has %s minor version lower then the application.' % what)
    if current[1] < loaded[1]:
        logger.warning('Loaded data has %s minor version higher then the application.' % what)
    if current[1] != loaded[1]:
        logger.warning('Possible problems with API compliance.')


def check_versions(version_dict):
    assert isinstance(version_dict, (dict, OrderedDict)), "File is corrupted, cannot read version data."
    assert 'version' in version_dict, "File is corrupted, cannot read version data."
    assert 'aquaduct_version' in version_dict, "File is corrupted, cannot read version data."
    check_version_compliance(aquaduct_version(), version_dict['aquaduct_version'], 'Aqua-Duct')
    check_version_compliance(version(), version_dict['version'], 'Valve')


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


def get_clustering_method(coptions,config):
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
        if abo.auto_barber:
            opts.update({'selection': str(abo.auto_barber)})
            if is_number(abo.auto_barber_mincut):
                opts.update({'mincut': float(abo.auto_barber_mincut)})
            else:
                opts.update({'mincut': None})
            if is_number(abo.auto_barber_maxcut):
                opts.update({'maxcut': float(abo.auto_barber_maxcut)})
            else:
                opts.update({'maxcut': None})
            opts.update({'mincut_level': bool(abo.auto_barber_mincut_level)})
            opts.update({'maxcut_level': bool(abo.auto_barber_maxcut_level)})
            opts.update({'tovdw': bool(abo.auto_barber_tovdw)})
        if 'auto_barber' in coptions._asdict():
            opts.update({'selection': str(coptions.auto_barber)})
        if 'auto_barber_maxcut' in coptions._asdict():
            if is_number(coptions.auto_barber_maxcut):
                opts.update({'maxcut': float(coptions.auto_barber_maxcut)})
            else:
                opts.update({'maxcut': None})
        if 'auto_barber_mincut' in coptions._asdict():
            if is_number(coptions.auto_barber_mincut):
                opts.update({'mincut': float(coptions.auto_barber_mincut)})
            else:
                opts.update({'mincut': None})
        if 'auto_barber_mincut_level' in coptions._asdict():
            opts.update({'mincut_level': bool(coptions.auto_barber_mincut_level)})
        if 'auto_barber_maxcut_level' in coptions._asdict():
            opts.update({'maxcut_level': bool(coptions.auto_barber_maxcut_level)})
        if 'auto_barber_tovdw' in coptions._asdict():
            opts.update({'tovdw': bool(coptions.auto_barber_tovdw)})

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
    clui.message('Goodby!')


def valve_load_config(filename, config):
    assert filename is not None, "No config file provided."
    assert os.path.isfile(filename), "Config file %s does not exist." % filename
    with clui.fbm('Load configuration file'):
        config.load_config(filename)


def valve_read_trajectory(top, traj, frames_window=None):
    with clui.fbm('Read trajectory'):
        return TrajectoryReader(top, traj, frames_window=frames_window)
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


def valve_exec_stage(stage, config, stage_run, reader=None, no_io=False, run_status=None,
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
            result = stage_run(config, options, reader=reader, **kwargs)
            run_status.update({stage: True})
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
                **kwargs):

    # create pool of workers - mapping function
    map_fun = map
    if optimal_threads.threads_count > 1:
        pool = mp.Pool(optimal_threads.threads_count)
        map_fun = pool.map

    with reader.get() as traj_reader:

        clui.message("Loop over frames - search of residues in object:")
        pbar = clui.pbar(traj_reader.number_of_frames)

        scope = traj_reader.parse_selection(options.scope)
        # scope will be used to derrive center of system
        center_of_system = np.array([0., 0., 0.])

        # create some containers
        res_ids_in_object_over_frames = {}
        all_res = None

        # the loop over frames
        for frame in traj_reader.iterate_over_frames():
            # center of system
            center_of_system += scope.center_of_mass()
            # current res selection
            res = traj_reader.parse_selection(options.object)
            # find matching residues:
            res_new = scope.containing_residues(res, convex_hull=options.scope_convexhull, map_fun=map_fun)
            # adds them to all_res
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
    if optimal_threads.threads_count > 1:
        pool.close()
        pool.join()
        del pool

    pbar.finish()
    center_of_system /= (frame + 1)
    logger.info('Center of system is %0.2f, %0.2f, %0.2f' % tuple(center_of_system))

    if all_res is None:
        raise ValueError("No traceable residues was found.")

    clui.message("Number of residues to trace: %d" % all_res.unique_resids_number())

    return {'all_res': all_res,
            'res_ids_in_object_over_frames': res_ids_in_object_over_frames,
            'center_of_system': center_of_system,
            'options': options}


################################################################################

# raw_paths
def stage_II_run(config, options,
                 reader=None,
                 all_res=None,
                 res_ids_in_object_over_frames=None,
                 **kwargs):
    if options.clear_in_object_info:
        clui.message('Clear data on residues in object over frames.')
        clui.message('This will be recalculated on demand.')
        res_ids_in_object_over_frames = {}

    with reader.get() as traj_reader:

        with clui.fbm("Init paths container"):
            paths = dict(((resid, GenericPaths(resid, min_pf=0, max_pf=traj_reader.number_of_frames-1)) for resid in
                          all_res.unique_resids(ikwid=True)))

        scope = traj_reader.parse_selection(options.scope)

        with clui.fbm("Rebuild treceable residues with current trajectory"):
            all_res = rebuild_selection(all_res, traj_reader)

        clui.message("Trajectory scan:")
        pbar = clui.pbar(traj_reader.number_of_frames)

        # create pool of workers - mapping function
        map_fun = map
        if optimal_threads.threads_count > 1:
            pool = mp.Pool(optimal_threads.threads_count)
            map_fun = pool.map

        for frame in traj_reader.iterate_over_frames():

            all_res_coords = list(all_res.center_of_mass_of_residues())  # this uses iterate over residues
            all_resids = [res.first_resid() for res in all_res.iterate_over_residues()]
            # check if is res are in scope
            is_res_in_scope = scope.contains_residues(all_res, convex_hull=options.scope_convexhull, map_fun=map_fun)

            # loop over coord, is  in scope, and resid
            for nr, (coord, isscope, resid) in enumerate(zip(all_res_coords, is_res_in_scope, all_resids)):
                # the point is that nr is not pointing to correct element in paths
                # now, nr is useless because paths is a dictionary, use resids instead
                assert paths[resid].id == resid, \
                    "Internal error. Paths IDs not synced with resids. \
                     Please send a bug report to the developer(s): %s" % __mail__
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

    if optimal_threads.threads_count > 1:
        pool.close()
        pool.join()
        del pool

    pbar.finish()

    clui.message("Number of paths: %d" % len(paths))

    return {'all_res': all_res, 'paths': paths, 'options': options}


################################################################################

class ABSphere(namedtuple('ABSphere', 'center radius')):
    pass


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
        with reader.get() as traj_reader:
            wtc = WhereToCut(spaths,traj_reader,
                             selection=options.auto_barber,
                             mincut=options.auto_barber_mincut,
                             mincut_level=options.auto_barber_mincut_level,
                             maxcut=options.auto_barber_maxcut,
                             maxcut_level=options.auto_barber_maxcut_level,
                             tovdw=options.auto_barber_mincut)
            # cut thyself!
            wtc.cut_thyself()

        clui.message("Auto Barber in action:")
        pbar = clui.pbar(len(paths))
        for p in paths.values():
            p.barber_with_spheres(wtc.spheres)
            pbar.next()
        pbar.finish()
        # now, it might be that some of paths are empty
        # paths = {k: v for k, v in paths.iteritems() if len(v.coords) > 0}
        paths = dict((k, v) for k, v in paths.iteritems() if
                     len(v.coords) > 0)  # more universal as dict comprehension may not work in <2.7
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
    operator_dict = {'>': operator.ge,
                     '=>': operator.gt,
                     '<=': operator.lt,
                     '<': operator.le}
    assert op in operator_dict.keys(), "Unsupported operator %s in threshold %s" % (op, rt)
    return lambda size_of_cluster: operator_dict[op](vl, size_of_cluster)

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

    def __init__(self,ths_def):

        self.thresholds = []
        if isinstance(ths_def,(str,unicode)):
            for thd in ths_def.split():
                self.thresholds.append(get_allow_size_function(thd))

    def __call__(self,size_of_cluster):
        for thd in self.thresholds:
            if not thd(size_of_cluster):
                return True
        return False

def potentially_recursive_clusterization(config=None,
                                         clusterization_name=None,
                                         inlets_object=None,
                                         spaths=None,
                                         traj_reader=None,
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
        clustering_function = get_clustering_method(cluster_options,config)
        # special case of barber!!!
        if cluster_options.method == 'barber':
            inlets_refs = inlets_object.get_inlets_references()
            wtc = WhereToCut([sp for sp in spaths if sp.id in inlets_refs],
                             traj_reader,
                             forceempty=True,
                             **clustering_function.method_kwargs)
            radii = [sphe.radius for sphe in wtc.spheres]
            inlets_object.add_radii(radii)
        # get skip_size function according to recursive_treshold
        #skip_size = get_skip_size_function(cluster_options.recursive_threshold)
        skip_size = SkipSizeFunction(cluster_options.recursive_threshold)
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
                                                    traj_reader=traj_reader,
                                                    deep=deep,
                                                    max_level=max_level)


# inlets_clusterization
def stage_IV_run(config, options,
                 spaths=None,
                 center_of_system=None,
                 reader=None,
                 **kwargs):
    coptions = config.get_cluster_options()
    rcoptions = config.get_recluster_options()
    soptions = config.get_smooth_options()

    max_level = int(options.max_level)
    assert max_level >= 0

    # new style clustering
    with clui.fbm("Create inlets"):
        # here we can check center of system
        inls = Inlets(spaths, center_of_system=center_of_system)
    clui.message("Number of inlets: %d" % inls.size)

    def noo():
        # returns number of outliers
        if 0 in inls.clusters_list:
            return inls.clusters.count(0)
        return 0

    if inls.size > 0:
        # ***** CLUSTERIZATION *****
        with clui.fbm("Performing clusterization", cont=False):
            with reader.get() as traj_reader:
                potentially_recursive_clusterization(config=config,
                                                     clusterization_name=config.cluster_name(),
                                                     inlets_object=inls,
                                                     spaths=spaths,
                                                     traj_reader=traj_reader,
                                                     message='clusterization',
                                                     max_level=max_level)
        # with log.fbm("Performing clusterization"):
        #    clustering_function = get_clustering_method(coptions)
        #    inls.perform_clustering(clustering_function)
        clui.message('Number of outliers: %d' % noo())
        # ***** OUTLIERS DETECTION *****
        if options.detect_outliers:
            with clui.fbm("Detecting outliers",cont=False):
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
                with reader.get() as traj_reader:
                    '''
                    potentially_recursive_clusterization(config=config,
                                                     clusterization_name=config.recluster_name(),
                                                     inlets_object=inls,
                                                     spaths=spaths,
                                                     traj_reader=traj_reader,
                                                     message='reclusterization',
                                                     max_level=max_level)
                    '''
                    clustering_function = get_clustering_method(rcoptions,config)
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

        clui.message('Clustering history:')
        clui.message(clui.print_simple_tree(inls.tree,prefix='').rstrip())

        with clui.fbm("Calculating cluster types"):
            ctypes = inls.spaths2ctypes(spaths)

        # now, there is something to do with ctypes!
        # we can create master paths!
        # but only if user wants this
        master_paths = {}
        master_paths_smooth = {}
        if options.create_master_paths:
            with clui.fbm("Creating master paths for cluster types", cont=False):
                smooth = get_smooth_method(soptions)
                ctypes_generic = [ct.generic for ct in ctypes]
                ctypes_generic_list = sorted(list(set(ctypes_generic)))

                pbar = clui.pbar(len(spaths) * 2)
                for nr, ct in enumerate(ctypes_generic_list):
                    logger.debug('CType %s (%d)' % (str(ct), nr))
                    sps = lind(spaths, what2what(ctypes_generic, [ct]))
                    logger.debug('CType %s (%d), number of spaths %d' % (str(ct), nr, len(sps)))
                    # print len(sps),ct
                    ctspc = CTypeSpathsCollection(spaths=sps, ctype=ct, pbar=pbar, threads=optimal_threads.threads_count)
                    master_paths.update({ct: ctspc.get_master_path(resid=nr)})
                    master_paths_smooth.update({ct: ctspc.get_master_path(resid=nr, smooth=smooth)})
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

    # TODO: Change it in such a way that it cooperates well with debug-file option.
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
        if len(t) > 1:  # traces.length_step_std requires at least 2 points
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
        # clui.message(sep())
        # clui.message('')
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

    '''
    # sorted by sizes:
    ctypes_generic_list = sorted(ctypes_generic_list, key=lambda ctyp: ctypes_size[ctypes_generic_list.index(ctyp)],
                                 reverse=True)
    '''

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
                scope = traj_reader.parse_selection(options.show_chull)
                frames_to_show = range2int(options.show_chull_frames)
                for frame in frames_to_show:
                    traj_reader.set_real_frame(frame)
                    chull = scope.get_convexhull_of_atom_positions()
                    spp.convexhull(chull, state=frame + 1)

    if options.show_object:
        with clui.fbm("Object shape"):
            with reader.get() as traj_reader:
                object_shape = traj_reader.parse_selection(options.show_object)
                frames_to_show = range2int(options.show_object_frames)
                for frame in frames_to_show:
                    traj_reader.set_real_frame(frame)
                    chull = object_shape.get_convexhull_of_atom_positions()
                    spp.convexhull(chull, name='object_shape', color=np.array([255, 153, 0]) / 255.,
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
version_nice
ValveConfig
valve_begin
valve_load_config
valve_read_trajectory
valve_exec_stage
stage_I_run
stage_II_run
stage_III_run
stage_IV_run
stage_V_run
stage_VI_run
valve_end
'''.split()
