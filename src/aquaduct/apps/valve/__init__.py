# -*- coding: utf-8 -*-

# This program is distributed in the hope that it will be useful,
# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2016-2018  Tomasz Magdziarz, Alicja Płuciennik, Michał Stolarczyk <info@aquaduct.pl>
# Copyright (C) 2019  Tomasz Magdziarz, Michał Banas <info@aquaduct.pl>
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

'''
This module is meant to read configuration file and init some global options.
It cannot import sandwich. Directly or indirectly.
'''

from aquaduct import logger, greetings as greetings_aquaduct, version_nice as aquaduct_version_nice

import configparser
import os
from collections import OrderedDict, namedtuple
from keyword import iskeyword

from aquaduct.apps.data import GCS, load_cric
from aquaduct.geom.cluster_available_methods import get_required_params, \
    AVAILABLE_METHODS as available_clustering_methods
from aquaduct.geom.smooth_available_methods import AVAILABLE_SMOOTHING_METHODS as available_smoothing_methods
from aquaduct.utils import clui
from aquaduct.utils.helpers import Auto


class ConfigSpecialNames(object):
    special_names_dict = {'none': None,
                          'null': None,
                          'true': True,
                          'false': False,
                          'yes': True,
                          'no': False,
                          'auto': Auto}

    def special_name(self, name):
        if isinstance(name, str):
            if name.lower() in self.special_names_dict:
                return self.special_names_dict[name.lower()]
        return name


class ValveConfig(ConfigSpecialNames):
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
        options_nt = namedtuple('Options', list(options.keys()))
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
        return 'scope scope_convexhull scope_everyframe scope_convexhull_inflate object'.split()

    @staticmethod
    def global_name():
        return 'global'

    @staticmethod
    def cluster_name():
        return 'clustering'

    @staticmethod
    def recluster_name():
        return 'reclustering'

    @staticmethod
    def recursive_clustering_name():
        return 'recursive_clustering'

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
                return 'inlets_clustering'
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
        # special recursive_clustering option
        if self.recursive_clustering_name() not in options:
            options.update({self.recursive_clustering_name(): None})
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

        config = configparser.RawConfigParser()

        def common(section):
            for setting in self.common_config_names():
                value = None
                if setting == 'execute':
                    value = 'runonce'
                elif setting == 'dump':
                    if GCS.netcdf:
                        value = '%d_%s_data.nc' % (snr + 1, section)
                    else:
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

        config.set(section, 'twoway', 'True')
        config.set(section, 'sps', 'True')

        config.set(section, 'sandwich', 'False')
        config.set(section, 'waterfall', 'False')
        config.set(section, 'max_frame', 'None')
        config.set(section, 'min_frame', 'None')
        config.set(section, 'step_frame', 'None')

        config.set(section, 'cache_dir', 'None')
        config.set(section, 'cache_mem', 'False')
        config.set(section, 'cache_type', 'full')

        ################
        snr = 0  # stage number
        # stage I
        # find traceable residues
        section = self.stage_names(snr)
        config.add_section(section)

        common(section)
        common_traj_data(section)
        config.set(section, 'scope_everyframe', 'False')
        config.set(section, 'scope_convexhull_inflate', 'None')

        config.set(section, 'add_passing', 'None')

        ################
        snr += 1
        # stage II
        # find raw paths
        section = self.stage_names(snr)
        config.add_section(section)

        common(section)
        common_traj_data(section)

        config.set(section, 'clear_in_object_info', 'False')
        config.set(section, 'discard_singletons', 1)
        config.set(section, 'discard_empty_paths', 'True')

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
        config.set(section, 'separate_barber', 'True')
        config.set(section, 'discard_empty_paths', 'True')
        config.set(section, 'sort_by_id', 'True')
        config.set(section, 'apply_smoothing', 'False')
        config.set(section, 'apply_soft_smoothing', 'True')
        config.set(section, 'remove_unused_parts', 'True')

        config.set(section, 'calculate_coo', 'False')

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
        config.set(section, 'exclude_passing_in_clustering', 'True')
        config.set(section, 'add_passing_to_clusters', 'None')
        config.set(section, 'renumber_clusters', 'False')
        config.set(section, 'join_clusters', 'None')
        config.set(section, 'master_paths_amount', 'None')
        config.set(section, 'separate_master', 'False')
        config.set(section, 'separate_master_all', 'True')
        config.set(section, 'inlets_center', 'cos')
        config.set(section, 'remove_inlets', 'None')
        config.set(section, 'clustering_order', 'old-school')

        ################
        # smooth
        section = self.smooth_name()
        config.add_section(section)
        config.set(section, 'method', 'window')

        ################
        # clustering
        section = self.cluster_name()
        config.add_section(section)
        config.set(section, 'method', 'barber')
        config.set(section, self.recursive_clustering_name(), self.cluster_name())
        config.set(section, self.recursive_threshold_name(), 'False')

        ################
        # reclustering
        section = self.recluster_name()
        config.add_section(section)
        config.set(section, 'method', 'dbscan')
        config.set(section, self.recursive_clustering_name(), 'False')
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
        config.set(section, 'scope_chull_inflate', 'None')
        config.set(section, 'object_chull', 'None')

        config.set(section, 'cluster_area', 'True')
        config.set(section, 'cluster_area_precision', '20')
        config.set(section, 'cluster_area_expand', '2')

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
        config.set(section, 'split_by_type', 'False')
        config.set(section, 'retain_all_types', 'False')

        # visualize spaths, all paths in one object
        config.set(section, 'all_paths_raw', 'False')
        config.set(section, 'all_paths_smooth', 'False')
        config.set(section, 'all_paths_split', 'False')  # split by in obj out
        config.set(section, 'all_paths_raw_io', 'False')
        config.set(section, 'all_paths_smooth_io', 'False')
        config.set(section, 'all_paths_amount', 'None')

        # visualize spaths, separate objects
        config.set(section, 'paths_raw', 'False')
        config.set(section, 'paths_smooth', 'False')
        config.set(section, 'paths_states', 'False')
        config.set(section, 'paths_raw_io', 'False')
        config.set(section, 'paths_smooth_io', 'False')

        config.set(section, 'ctypes_raw', 'False')
        config.set(section, 'ctypes_smooth', 'False')
        config.set(section, 'ctypes_amount', 'None')

        # visualize clusters
        config.set(section, 'inlets_clusters', 'False')
        config.set(section, 'inlets_clusters_amount', 'None')

        config.set(section, 'cluster_area', 'False')
        config.set(section, 'cluster_area_precision', '20')
        config.set(section, 'cluster_area_expand', '2')

        # show protein
        config.set(section, 'show_molecule', 'None')
        config.set(section, 'show_molecule_frames', '0')
        config.set(section, 'show_scope_chull', 'None')
        config.set(section, 'show_scope_chull_inflate', 'None')
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
            out = ['Possible clustering methods:']
            # get default clustering method
            defmet = self.config.get(section, 'method')
            for method in available_clustering_methods:
                if method == defmet:
                    continue
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
        if section == self.smooth_name():
            out = ['Possible smoothing methods:']
            # get default clustering method
            defmet = self.config.get(section, 'method')
            for method in available_smoothing_methods:
                if method == defmet:
                    continue
                out.append('method = %s' % method)
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
            for key, value in opts._asdict().items():  # loop over options
                if key in skip_list:
                    continue
                # comment scope etc. in stage II if dump_template
                if dump_template and name == self.stage_names(1):
                    if key in self.common_traj_data_config_names():
                        key = '#' + key
                output.append('%s = %s' % (key, value2str(value)))
            if not concise:
                output.append('')

        # is something missing? another loop over all additional sections
        for miss in self.config.sections():
            if miss in names:
                continue  # skip if it was already dumped in the loop above
            output.append('[%s]' % miss)
            for key in self.config.options(miss):  # loop over options in section miss
                if key in skip_list:
                    continue
                value = self.config.get(miss, key)
                output.append('%s = %s' % (key, value2str(value)))
            if not concise:
                output.append('')
        while not len(output[-1].strip()):
            output.pop()
        return output


def valve_load_config(filename, config):
    assert filename is not None, "No config file provided."
    assert os.path.isfile(filename), "Config file %s does not exist." % filename
    with clui.fbm('Load configuration file'):
        config.load_config(filename)
        goptions = config.get_global_options()

        # coords cache type
        if goptions.cache_type not in ['full','simple']:
            logger.warning("Unknown cache type '%s', using 'full' by default" % goptions.cache_type)
        else:
            logger.info("Using '%s' cache type." % goptions.cache_type)
        GCS.cachetype = goptions.cache_type
        # cache dir & netcdf
        GCS.cachedir = goptions.cache_dir
        GCS.cachemem = goptions.cache_mem
        # GCS.netcdf = args.use_netcdf or args.use_netcdf4
        # GCS.nc4 = args.use_netcdf4
        load_cric()

        # single precision storage
        if goptions.sps:
            logger.info('Single precision data storage activated.')
            from aquaduct.utils.maths import defaults
            from numpy import float32, int32

            defaults.float_default = float32
            defaults.int_default = int32


def sep():
    return clui.gsep(sep='-', times=48)


def asep():
    return clui.gsep(sep='=', times=72)


def valve_begin():
    clui.message(greetings_aquaduct())  # nice greetings
    clui.message('Aqua-Duct version %s' % aquaduct_version_nice())
    # clui.message('Valve driver version %s' % version_nice())
    clui.message(sep())


def valve_end():
    clui.message(sep())
    clui.message('Let the Valve be always open!')
    clui.message('Goodbye!')
