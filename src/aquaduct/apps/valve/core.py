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

import ConfigParser
import os
import gc
from collections import namedtuple, OrderedDict
from functools import partial
from itertools import izip_longest, izip, chain
from keyword import iskeyword
from multiprocessing import Pool, Queue, Process

from scipy.spatial.distance import cdist

from aquaduct import greetings as greetings_aquaduct
from aquaduct import logger
from aquaduct import version_nice as aquaduct_version_nice
from aquaduct.apps.data import GCS, CRIC, save_cric
from aquaduct.apps.valve.data import get_vda_reader
from aquaduct.geom.cluster import AVAILABLE_METHODS as available_clusterization_methods
from aquaduct.geom.cluster import get_required_params
from aquaduct.geom.master import CTypeSpathsCollection
from aquaduct.traj.barber import WhereToCut, barber_paths
from aquaduct.traj.inlets import InletClusterGenericType
from aquaduct.traj.inlets import Inlets
from aquaduct.traj.paths import GenericPaths, yield_single_paths, SinglePath
from aquaduct.traj.paths import yield_generic_paths
from aquaduct.traj.sandwich import ResidueSelection, Reader
from aquaduct.utils import clui
from aquaduct.utils.clui import roman
from aquaduct.utils.helpers import iterate_or_die, fractionof
from aquaduct.utils.helpers import range2int, Auto, what2what, lind, is_number, robust_and, robust_or
from aquaduct.utils.multip import optimal_threads
from aquaduct.utils.sets import intersection_full

from aquaduct.apps.valve.worker import stage_I_worker_q, stage_II_worker_q, stage_II_worker_q_twoways, \
    assign_nonsandwiched_paths, assign_sandwiched_paths
from aquaduct.apps.valve.helpers import *
from aquaduct.apps.valve.spath import *
from aquaduct.apps.valve.clusters import *

__mail__ = 'info@aquaduct.pl'
__version__ = aquaduct_version_nice()


###############################################################################
# configuration file helpers

class ConfigSpecialNames(object):
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
        return 'scope scope_convexhull scope_everyframe scope_convexhull_inflate object'.split()

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
        config.set(section, 'scope_chull_inflate', 'None')
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
            out = ['Possible clusterization methods:']
            # get default clusterization method
            defmet = self.config.get(section, 'method')
            for method in available_clusterization_methods:
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


def valve_exec_stage(stage, config, stage_run, no_io=False, run_status=None, force_save=None,
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
                    vda = get_vda_reader(options.dump, mode='w')
                    vda.dump(**result)
                # save_stage_dump(options.dump, **result)
        elif options.execute in ['skip'] or (options.execute in ['runonce'] and can_be_loaded):
            if not no_io:
                ###########
                # L O A D #
                ###########
                if options.dump:
                    with clui.fbm('Loading data dump from %s file' % options.dump):
                        vda = get_vda_reader(options.dump, mode='r')
                        result = vda.load()
                        # result = load_stage_dump(options.dump, reader=reader)
                    if force_save:
                        ######################
                        # F O R C E  S A V E #
                        ######################
                        if force_save in ['dump', 'nc']:
                            fs_fname = os.path.extsep.join([os.path.splitext(options.dump)[0], force_save])
                        else:
                            fs_fname = options.dump
                        with clui.fbm('Saving (forced) data dump in %s file' % fs_fname):
                            vda = get_vda_reader(fs_fname, mode='w')
                            vda.dump(**result)
        else:
            raise NotImplementedError('exec mode %s not implemented' % options.execute)
        # remove options stuff
        # GC!
        gc.collect()
        if not no_io:
            if result is not None:
                return dict(((key, val) for key, val in result.iteritems() if 'options' not in key))


################################################################################
# stages run


# traceable_residues
def stage_I_run(config, options,
                **kwargs):
    Reader.reset()

    clui.message("Loop over frames - search of residues in object:")
    pbar = clui.pbar(Reader.number_of_frames())

    # create some containers
    number_frame_rid_in_object = []
    # res_ids_in_scope_over_frames = {}  # not used
    all_res = None
    # scope will be used to derrive center of system
    center_of_system = np.array([0., 0., 0.])

    # prepare queues
    pbar_queue = Queue()
    results_queue = Queue()
    input_queue = Queue()

    # prepare and start pool of workers
    pool = [Process(target=stage_I_worker_q, args=(input_queue, results_queue, pbar_queue)) for dummy in
            xrange(optimal_threads.threads_count)]
    [p.start() for p in pool]

    if options.scope_convexhull_inflate:
        logger.info("Inflate convex hull by %0.1f" % float(options.scope_convexhull_inflate))
    # feed input_queue with data
    for results_count, (number, traj_reader) in enumerate(Reader.iterate(number=True)):
        input_queue.put((number, traj_reader,
                         options.scope_everyframe, options.scope, options.scope_convexhull, options.scope_convexhull_inflate,
                         options.object, max(1, Reader.number_of_frames() / 500)))

    # display progress
    progress = 0
    progress_target = Reader.number_of_frames()
    for p in iter(pbar_queue.get, None):
        pbar.next(p)
        progress += p
        if progress == progress_target:
            break
    # [stop workers]
    [input_queue.put(None) for p in pool]
    pbar.finish()

    # collect results
    pbar = clui.pbar((results_count + 1) * 2 + 1, 'Collecting results from layers:')
    results = {}
    for nr, result in enumerate(iter(results_queue.get, None)):
        results.update(result)
        pbar.next()
        if nr == results_count:
            break
    for key in sorted(results.keys()):
        _all_res, _frame_rid_in_object, _center_of_system = results.pop(key)
        center_of_system += _center_of_system
        if all_res:
            all_res.add(_all_res)
            all_res.uniquify()
        else:
            all_res = _all_res
        number_frame_rid_in_object.append(_frame_rid_in_object)
        pbar.next()
    center_of_system /= (Reader.number_of_frames())
    logger.info('Center of system is %0.2f, %0.2f, %0.2f' % tuple(center_of_system))
    pbar.next()
    # join pool
    [p.join(1) for p in pool]
    pbar.finish()

    if all_res is None or all_res.len() == 0:
        raise ValueError("No traceable residues were found.")

    if options.add_passing:
        with clui.fbm("Add possible passing paths mols"):
            for traj_reader in Reader.iterate():
                traj_reader = traj_reader.open()
                passing = traj_reader.parse_selection(options.add_passing).residues()
                all_res.add(passing)
                all_res.uniquify()

    unsandwitchize = not Reader.sandwich_mode and len(number_frame_rid_in_object) > 1
    if unsandwitchize:
        with clui.fbm("Unsandwitchize traced residues"):
            # all_res, each layer should comprise of the same res
            all_res_ids = sorted(set([i[-1] for i in all_res.ids()]))
            all_res = ResidueSelection({0: all_res_ids})
            # do similar operation with number_frame_rid_in_object
            while (len(number_frame_rid_in_object)) > 1:
                number_frame_rid_in_object[0].extend(number_frame_rid_in_object.pop(1))

    clui.message("Number of residues to trace: %d" % all_res.len())

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
    # disable real cache of ort
    Reader.reset()

    ####################################################################################################################

    # clear in object info if required
    if options.clear_in_object_info:
        clui.message('Clear data on residues in object over frames.')
        clui.message('This will be recalculated on demand.')
        number_frame_rid_in_object = None

    is_number_frame_rid_in_object = bool(number_frame_rid_in_object)

    number_of_frames = Reader.number_of_frames(onelayer=True)

    with clui.pbar(Reader.number_of_frames(), mess="Trajectory scan:") as pbar:

        # prepare queues
        pbar_queue = Queue()
        results_queue = Queue()
        input_queue = Queue()

        allow_twoway = False
        if config.get_global_options().twoway:
            sI = config.get_stage_options(0)
            sII = config.get_stage_options(1)
            sIII = config.get_stage_options(2)
            # TODO: enable twoway in sandwich mode
            allow_twoway = (not Reader.sandwich_mode) and (sI.scope == sII.scope) and (sI.scope_convexhull == sII.scope_convexhull) and (sI.object == sII.object) and (sIII.allow_passing_paths == False)

        # prepare and start pool of workers
        if allow_twoway:
            logger.info('Twoway trajectory scan enabled.')
            pool = [Process(target=stage_II_worker_q_twoways, args=(input_queue, results_queue, pbar_queue)) for dummy in
                    xrange(optimal_threads.threads_count)]
        else:
            pool = [Process(target=stage_II_worker_q, args=(input_queue, results_queue, pbar_queue)) for dummy in
                    xrange(optimal_threads.threads_count)]
        [p.start() for p in pool]

        if options.scope_convexhull_inflate:
            logger.info("Inflate convex hull by %0.1f" % float(options.scope_convexhull_inflate))
        # feed input_queue with data
        if Reader.sandwich_mode:
            for results_count, (frame_rid_in_object, (number, traj_reader)) in enumerate(
                    izip(iterate_or_die(number_frame_rid_in_object,
                                        times=Reader.number_of_layers()),
                         Reader.iterate(number=True))):
                all_res_layer = all_res.layer(number)

                input_queue.put((number, traj_reader,
                                 options.scope_everyframe, options.scope, options.scope_convexhull, options.scope_convexhull_inflate,
                                 options.object, all_res_layer, frame_rid_in_object, is_number_frame_rid_in_object,
                                 max(1, optimal_threads.threads_count)))
        else:
            all_res_ids = sorted(set([i[-1] for i in all_res.ids()]))
            seek = 0
            frame_rid_in_object = None
            for results_count, (number, traj_reader) in enumerate(Reader.iterate(number=True)):
                all_res_layer = ResidueSelection({number: all_res_ids})
                if is_number_frame_rid_in_object:
                    frame_rid_in_object = number_frame_rid_in_object[0][seek:seek + traj_reader.window.len()]
                    seek += traj_reader.window.len()
                input_queue.put((number, traj_reader,
                                 options.scope_everyframe, options.scope, options.scope_convexhull, options.scope_convexhull_inflate,
                                 options.object, all_res_layer, frame_rid_in_object, is_number_frame_rid_in_object,
                                 max(1, optimal_threads.threads_count)))
            all_res_names = list(all_res_layer.names())

        # display progress
        progress = 0
        progress_target = Reader.number_of_frames()
        # progress_target = len(all_res)
        for p in iter(pbar_queue.get, None):
            pbar.next(p)
            progress += p
            if progress == progress_target:
                break
            if progress % (max(1, optimal_threads.threads_count)**2 * 1000) == 0:
                gc.collect()

        # [stop workers]
        [input_queue.put(None) for p in pool]

    # collect results
    pbar = clui.pbar((results_count + 1) + 1, 'Collecting results from layers:')
    results = {}
    for nr, result in enumerate(iter(results_queue.get, None)):
        for rk in result.keys():
            if rk not in results:
                results.update({rk: result[rk]})
            else:
                raise ValueError('Wrong results format; please send bug report.')
                results.update({rk: results[rk] + result[rk]})
        pbar.next()
        if nr == results_count:
            break
    [p.join(1) for p in pool]
    pbar.next()
    pbar.finish()

    # now, results holds 012 matrices, make paths out of it
    # TODO: Following could be easily written in parallel.

    unsandwitchize = not Reader.sandwich_mode and len(results) > 1
    if not unsandwitchize:
        with clui.pbar(len(all_res), 'Creating raw paths:') as pbar:
            paths = []
            for number in sorted(results.keys()):
                all_res_layer = all_res.layer(number)

                paths_this_layer = (GenericPaths(resid,
                                                 name_of_res=resname,
                                                 min_pf=0, max_pf=number_of_frames - 1)
                                    for resid, resname in izip(all_res_layer.ids(),
                                                               all_res_layer.names()))

                pool = Pool(processes=optimal_threads.threads_count)
                r = pool.map_async(
                    assign_nonsandwiched_paths(pbar),
                    izip(paths_this_layer, results_n(results[number]).T),
                    callback=paths.extend)
                r.wait()

    elif unsandwitchize:
        # make coherent paths
        # paths names, paths
        max_pf = Reader.number_of_frames() - 1
        # frames_offset = np.cumsum([0] + frames_offset).tolist()[:len(numbers)]

        # pbar = clui.pbar(len(results[numbers[0]]), 'Sandwich deconvolution:')
        with clui.pbar(len(all_res_ids), 'Creating raw paths (sandwich deconvolution):') as pbar:
            paths = []

            pool_func = assign_sandwiched_paths(all_res_ids, all_res_names, max_pf, results, pbar)

            pool = Pool(processes=optimal_threads.threads_count)
            r = pool.map_async(
                pool_func, xrange(len(all_res_ids)),
                callback=paths.extend)
            r.wait()

    # rm tmp files
    for rn in results.itervalues():
        if not isinstance(rn, np.ndarray):
            os.unlink(rn[0])

    if options.discard_empty_paths:
        with clui.pbar(maxval=len(paths), mess="Discard residues with empty paths:") as pbar:
            paths = [pat for pat in paths if len(pat.frames) > 0 if pbar.next() is None]

    clui.message("Number of paths: %d" % len(paths))

    return {'all_res': all_res, 'paths': paths, 'options': options._asdict()}


################################################################################


# separate_paths
def stage_III_run(config, options,
                  paths=None,
                  **kwargs):
    soptions = config.get_smooth_options()

    Reader.reset()

    if options.allow_passing_paths:
        logger.warning("Passing paths is a experimental feature. Please, analyze results with care.")

    ######################################################################

    if options.discard_empty_paths:
        with clui.pbar(maxval=len(paths), mess="Discard residues with empty paths:") as pbar:
            paths = [pat for pat in paths if len(pat.frames) > 0 if pbar.next() is None]

    ######################################################################

    with clui.pbar(len(paths), "Create separate paths:") as pbar:
        Reader.reset()
        pool = Pool(processes=optimal_threads.threads_count)
        ysp = partial(yield_single_paths, progress=True, passing=options.allow_passing_paths)
        n = max(1, optimal_threads.threads_count)
        spaths = []
        nr_all = 0
        for sps_nrs in pool.imap_unordered(ysp, (paths[i:i + n] for i in xrange(0, len(paths), n))):
            nr = 0  # if no spaths were returned
            for sp, nr in sps_nrs:
                spaths.append(sp)
                pbar.update(nr_all + nr + 1)
            nr_all += nr + 1
        pool.close()
        pool.join()
        gc.collect()

    clui.message("Created %d separate paths out of %d raw paths" %
                 (len(spaths), len(paths)))

    ######################################################################

    with clui.pbar(len(spaths) + len(paths), "Removing unused parts of paths:") as pbar:
        paths = yield_generic_paths(spaths, progress=pbar)

    ######################################################################

    if options.discard_short_paths or options.discard_short_object:
        # prepare options
        if is_number(options.discard_short_paths):
            short_paths = int(options.discard_short_paths)
        else:
            short_paths = None
        if is_number(options.discard_short_object):
            short_object = float(options.discard_short_object)
        else:
            short_object = None
        # check logic
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
        # do calculations
        if short_paths is not None or short_object is not None:
            with clui.pbar(len(spaths), discard_message) as pbar:
                spaths_nr = len(spaths)
                Reader.reset()
                pool = Pool(processes=optimal_threads.threads_count)
                dse = partial(discard_short_etc, short_paths=short_paths, short_object=short_object,
                              short_logic=short_logic)
                n = max(1, optimal_threads.threads_count)
                spaths_new = pool.imap_unordered(dse, (spaths[i:i + n] for i in xrange(0, len(spaths), n)))
                # spaths_new = imap(dse, (spaths[i:i + n] for i in xrange(0, len(spaths), n)))
                # CRIC AWARE MP!
                if short_object is not None:
                    Reader.reset()
                    spaths = list(chain.from_iterable((sps for nr, sps, cric in spaths_new if
                                                       (pbar.next(step=nr) is None) and (
                                                               CRIC.update_cric(cric) is None))))
                    save_cric()
                else:
                    spaths = list(chain.from_iterable((sps for nr, sps in spaths_new if pbar.next(step=nr) is None)))

                pool.close()
                pool.join()
                gc.collect()
            # del spaths
            spaths_nr_new = len(spaths)
            if spaths_nr == spaths_nr_new:
                clui.message("No paths were discarded.")
            else:
                clui.message("%d paths were discarded." % (spaths_nr - spaths_nr_new))
                with clui.pbar(len(spaths) + len(paths), "Removing (again) unused parts of paths:") as pbar:
                    paths = yield_generic_paths(spaths, progress=pbar)
        else:
            clui.message("No paths were discarded - no values were set.")

    ######################################################################

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

        if len(wtc.spheres):
            with clui.pbar(maxval=len(paths), mess="AutoBarber in action:") as pbar:
                Reader.reset()
                pool = Pool(processes=optimal_threads.threads_count)
                bp = partial(barber_paths, spheres=wtc.spheres)
                n = max(1, optimal_threads.threads_count)
                paths_new = pool.imap_unordered(bp, (paths[i:i + n] for i in xrange(0, len(paths), n)))
                # paths_new = map(bp, (paths[i:i + n] for i in xrange(0, len(paths), n)))
                paths_ = list()
                for paths_new_list in paths_new:
                    CRIC.update_cric(paths_new_list.pop(-1))
                    paths_.extend(paths_new_list)
                    pbar.next(step=len(paths_new_list))
                save_cric()
                # now, it might be that some of paths are empty
                paths = []
                while len(paths_):
                    pat = paths_.pop(0)
                    if len(pat.frames) > 0:
                        paths.append(pat)
                del paths_,paths_new,paths_new_list
                pool.close()
                pool.join()
                gc.collect()
        else:
            clui.message('AutoBarber procedure skip, no spheres detected.')

    ######################################################################
    # following procedures are run only if autobarber
    ######################################################################

    if options.auto_barber:

        clui.message("Recreate separate paths:")
        pbar = clui.pbar(len(paths))
        Reader.reset()
        pool = Pool(processes=optimal_threads.threads_count)
        # ysp = partial(yield_single_paths, progress=True, passing=options.allow_passing_paths)
        n = max(1, len(paths) / optimal_threads.threads_count / 3)
        spaths = []
        nr_all = 0
        for sps_nrs in pool.imap_unordered(ysp, (paths[i:i + n] for i in xrange(0, len(paths), n))):
            for sp, nr in sps_nrs:
                spaths.append(sp)
                pbar.update(nr_all + nr + 1)
            nr_all += nr + 1
        pool.close()
        pool.join()
        gc.collect()
        pbar.finish()

    clui.message("(Re)Created %d separate paths out of %d raw paths" %
                 (len(spaths), len(paths)))


    ######################################################################

    if options.auto_barber:

        if options.discard_short_paths or options.discard_short_object:
            if short_paths is not None or short_object is not None:
                with clui.pbar(len(spaths), discard_message) as pbar:
                    spaths_nr = len(spaths)
                    Reader.reset()
                    pool = Pool(processes=optimal_threads.threads_count)
                    #dse = partial(discard_short_etc, short_paths=short_paths, short_object=short_object,
                    #              short_logic=short_logic)
                    n = max(1, optimal_threads.threads_count)
                    spaths_new = pool.imap_unordered(dse, (spaths[i:i + n] for i in xrange(0, len(spaths), n)))
                    # CRIC AWARE MP!
                    if short_object is not None:
                        Reader.reset()
                        spaths = list(chain.from_iterable((sps for nr, sps, cric in spaths_new if
                                                           (pbar.next(step=nr) is None) and (
                                                                   CRIC.update_cric(cric) is None))))
                    else:
                        spaths = list(
                            chain.from_iterable((sps for nr, sps in spaths_new if pbar.next(step=nr) is None)))
                    save_cric()
                    pool.close()
                    pool.join()
                    gc.collect()
                # del spathsqq
                spaths_nr_new = len(spaths)
                if spaths_nr == spaths_nr_new:
                    clui.message("No paths were discarded.")
                else:
                    clui.message("%d paths were discarded." % (spaths_nr - spaths_nr_new))
                    with clui.pbar(len(spaths) + len(paths), "Removing (again) unused parts of paths:") as pbar:
                        paths = yield_generic_paths(spaths, progress=pbar)
            else:
                clui.message("No paths were discarded - no values were set.")

    ######################################################################

    if options.sort_by_id:
        with clui.fbm("Sort separate paths by resid"):
            spaths = sorted(spaths, key=lambda sp: (sp.id.id, sp.id.nr))

    ######################################################################

    # apply smoothing?
    # it is no longer necessary
    if options.apply_smoothing:
        logger.warning(
            "Hard smoothing is not available in the current version but may be available in the future. Stay tuned!")
    if options.apply_soft_smoothing:
        logger.warning(
            "Soft smoothing option is not available any more. Soft smoothing is enabled by default if --cache-dir or --cache-mem options are used.")
    clui.message("Number of paths: %d" % len(paths))
    clui.message("Number of spaths: %d" % len(spaths))

    return {'paths': paths, 'spaths': spaths, 'options': options._asdict(), 'soptions': soptions._asdict()}


################################################################################


# inlets_clusterization
def stage_IV_run(config, options,
                 spaths=None,
                 center_of_system=None,
                 **kwargs):
    # enable real cache of ort
    Reader.reset()

    coptions = config.get_cluster_options()
    rcoptions = config.get_recluster_options()
    soptions = config.get_smooth_options()

    max_level = int(options.max_level)
    assert max_level >= 0

    # new style clustering
    # with clui.fbm("Create inlets"):
    # here we can check center of system
    pbar = clui.SimpleProgressBar(maxval=len(spaths), mess="Create inlets")
    inls = Inlets(spaths, center_of_system=center_of_system, passing=not options.exclude_passing_in_clusterization,
                  pbar=pbar)
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
        gc.collect()
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
        gc.collect()
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
        gc.collect()
        # ***** SINGLETONS REMOVAL *****
        if options.singletons_outliers:
            with clui.fbm("Removing clusters of size %d" % int(options.singletons_outliers)):
                inls.small_clusters_to_outliers(int(options.singletons_outliers))
            clui.message('Number of clusters detected so far: %d' % len(inls.clusters_list))
            clui.message('Number of outliers: %d' % noo())

        # TODO: Move it after master paths!
        gc.collect()
        # ***** ADD PASSING PATHS TO CLUSTERS *****
        if options.exclude_passing_in_clusterization and options.add_passing_to_clusters:
            with clui.fbm("Adding passing paths inlets to clusters", cont=False):
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
                        if cluster == 0:
                            continue
                        clui.message("Current cluster: %d." % cluster)
                        # sps = inls.lim2clusters(cluster).limspaths2(spaths_single)
                        # chull = inls.lim2clusters(cluster).get_chull()
                        wtc = WhereToCut(inlets=inls.lim2clusters(cluster), **ab_options)
                        wtc.cut_thyself()
                        pbar = clui.SimpleProgressBar(len(passing_inlets_ids),
                                                      "Loop over available passing paths inlets:")
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

        gc.collect()

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
                    pbar = clui.pbar(len(spaths) * 2, mess='Building coords cache')
                    # TODO: do it in parallel
                    [sp.get_coords(smooth=None) for sp in spaths if
                     pbar.next() is None and not isinstance(sp, PassingPath)]
                    [sp.get_coords(smooth=smooth) for sp in spaths if
                     pbar.next() is None and not isinstance(sp, PassingPath)]
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
                        sps = [sp for sp in sps if not isinstance(sp, PassingPath)]  # no PassingPaths!
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
        gc.collect()
    else:
        clui.message("No inlets found. Clusterization skipped.")
        # make empty results
        ctypes = inls.spaths2ctypes(spaths)
        master_paths = {}
        master_paths_smooth = {}

    ################################################################################

    return {'inls': inls,
            'ctypes': ctypes,
            'master_paths': master_paths,
            'master_paths_smooth': master_paths_smooth}


################################################################################


################################################################################

# analysis
def stage_V_run(config, options,
                spaths=None,
                paths=None,
                inls=None,
                ctypes=None,
                reader=None,
                **kwargs):
    # enable real cache of ort
    Reader.reset()

    # file handle?
    head_nr = True
    line_nr = head_nr
    pa = PrintAnalysis(options.save, line_nr=line_nr)

    if options.save:
        clui.message('Using user provided file (%s), and' % options.save)
        clui.message('for histograms data file (%s).' % (options.save + '.csv'))
        pbar = True
    else:
        clui.message('Using standard output.')
        clui.message(sep())
        clui.message('')
        pbar = False

    ############
    pa.sep()
    pa('Aqua-Duct analysis')
    pa(clui.get_str_timestamp())

    ############
    # format for path ID
    max_ID_len = 0
    for sp in spaths:
        n = len(str(sp.id))
        if n > max_ID_len:
            max_ID_len = n
    spath_id_header.format = '%%%ds' % (max_ID_len + 2)

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

    if pbar:
        # calculate number of tnspt
        nr_tnspt = [nr for nr, dummy in enumerate(iter_over_tnspt())][-1] + 1
        nr_f = (1 if len(traced_names) == 1 else 2) * (1 if len(spaths_types) == 1 else 2)

        pbar = clui.pbar(maxval=nr_tnspt * 2 + nr_f * 2 * len(spaths) + len(spaths), mess="Calculating stats:")

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
        if pbar:
            pbar.next()

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
            if pbar:
                pbar.heartbeat()
        pa.tend(header_line)

        header_line, line_template = get_header_line_and_line_template(clusters_stats_len_header(), head_nr=head_nr)
        pa("Clusters statistics (of paths%s) mean lengths of transfers" % message)
        pa.thead(header_line)
        for nr, cl in enumerate(inls.clusters_list):
            sp_ct_lim = ((sp, ct) for sp, ct in zip(spaths, ctypes) if
                         cl in ct.clusters and isinstance(sp, sptype) and sp.id.name in tname)
            pa(make_line(line_template, clusters_stats_len(cl, sp_ct_lim)), nr=nr)
            if pbar:
                pbar.heartbeat()
        pa.tend(header_line)

        header_line, line_template = get_header_line_and_line_template(clusters_stats_steps_header(), head_nr=head_nr)
        pa("Clusters statistics (of paths%s) mean frames numbers of transfers" % message)
        pa.thead(header_line)
        for nr, cl in enumerate(inls.clusters_list):
            sp_ct_lim = ((sp, ct) for sp, ct in zip(spaths, ctypes) if
                         cl in ct.clusters and isinstance(sp, sptype) and sp.id.name in tname)
            pa(make_line(line_template, clusters_stats_steps(cl, sp_ct_lim)), nr=nr)
            if pbar:
                pbar.heartbeat()
        pa.tend(header_line)
        if pbar:
            pbar.next()

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
            if pbar:
                pbar.next(len(sps))
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
            if pbar:
                pbar.next(len(sps))
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
        if pbar:
            pbar.next()

    pa.tend(header_line)

    if pbar:
        pbar.finish()

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
                     mess='Calculating histograms:')
    # loop over paths and ctypes
    for sp, ct in zip(spaths, ctypes):
        # loop over columns
        # traced names, paths types, clusters cluster types, part of paths, column name
        for tname, sptype, c_ct, part, col_name in iter_over_all():
            # check if column fits to the requirements
            if not sp.id.name in tname:
                continue
            if not isinstance(sp, sptype):
                continue
            # c or ct
            it_is_ct = False
            if isinstance(c_ct[0], InletClusterGenericType):
                it_is_ct = True
                if not ct.generic in c_ct:
                    continue
            else:
                if len(intersection_full(ct.generic.clusters, c_ct)) == 0:
                    continue
            # this is more than that...
            # if c_ct is not InletClusterGenericType then:
            # if 'in' in part only incoming paths in ct are used
            # if 'out' in part only outgoing paths in ct are used
            col_index = header.index(col_name)
            if 'walk' in part:
                h[sp.paths_cont, col_index] += 1
            if isinstance(sp, PassingPath):
                continue
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


    if options.scope_chull_inflate:
        logger.info("Inflate convex hull by %0.1f" % float(options.scope_chull_inflate))

    # scope/object size?
    if options.calculate_scope_object_size:
        scope_size = []
        object_size = []
        # now, the problem is in the scope and object definition.
        pbar = clui.pbar(maxval=Reader.number_of_frames(), mess='Calculating scope and object sizes')

        for number, traj_reader in Reader.iterate(number=True):
            traj_reader = traj_reader.open()
            if Reader.sandwich_mode:
                header += map(lambda s: '%s_%d' % (s, number),
                              ['scope_area', 'scope_volume', 'object_area', 'object_volume'])
                fmt += ['%0.3f', '%0.2f'] * 2
            elif not len(scope_size):
                header += map(lambda s: '%s' % s,
                              ['scope_area', 'scope_volume', 'object_area', 'object_volume'])
                fmt += ['%0.3f', '%0.2f'] * 2
            if Reader.sandwich_mode or not len(scope_size):
                scope_size.append([])
                object_size.append([])
            for frame in traj_reader.iterate_over_frames():
                scope = traj_reader.parse_selection(options.scope_chull)
                ch = scope.chull(inflate=options.scope_chull_inflate)
                scope_size[-1].append((ch.area, ch.volume))
                res = traj_reader.parse_selection(options.object_chull)
                ch = res.chull()
                object_size[-1].append((ch.area, ch.volume))
                pbar.next()
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


def stage_VI_run(config, options,
                 spaths=None,
                 inls=None,
                 ctypes=None,
                 master_paths=None,
                 master_paths_smooth=None,
                 **kwargs):
    # enable real cache of ort
    Reader.reset()

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
            for nr, traj_reader in enumerate(Reader.iterate(threads=False)):
                traj_reader = traj_reader.open()
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
            if options.show_scope_chull_inflate:
                logger.info("Inflate convex hull by %0.1f" % float(options.show_scope_chull_inflate))
            for nr, traj_reader in enumerate(Reader.iterate(threads=False)):
                traj_reader = traj_reader.open()
                frames_to_show = range2int(options.show_scope_chull_frames)
                for frame in frames_to_show:
                    traj_reader.set_frame(frame)
                    scope = traj_reader.parse_selection(options.show_scope_chull)
                    chull = scope.chull(inflate=options.show_scope_chull_inflate)
                    spp.convexhull(chull, name='scope_shape%d' % nr, state=frame + 1)

    if options.show_object_chull:
        with clui.fbm("Object shape"):
            for nr, traj_reader in enumerate(Reader.iterate(threads=False)):
                traj_reader = traj_reader.open()
                frames_to_show = range2int(options.show_object_chull_frames)
                for frame in frames_to_show:
                    traj_reader.set_frame(frame)
                    object_shape = traj_reader.parse_selection(options.show_object_chull)
                    chull = object_shape.chull()
                    spp.convexhull(chull, name='object_shape%d' % nr, color=np.array([255, 153, 0]) / 255.,
                                   state=frame + 1)  # orange


    def make_fracion(frac,size):
        if frac is not None:
            frac = float(frac)
            if frac > 1:
                frac = frac / size
                if frac >= 1:
                    frac = None
        return frac

    fof = lambda sp: np.array(list(fractionof(sp, f=make_fracion(options.inlets_clusters_amount,len(sp)))))

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
                spp.scatter(fof(ics), color=cmap(c), name="cluster_%s" % c_name)
                if False:  # TODO: This does not work any more in that way. Rewrite it or remove it
                    radii = inls.lim2clusters(c).radii
                    if len(radii) > 0:
                        spp.scatter(ics, color=cmap(c), radius=radii, name="cluster_radii_%s" % c_name)

    fof = lambda sp: np.array(list(fractionof(sp, f=make_fracion(options.ctypes_amount,len(sp)))))

    if options.ctypes_raw:
        with clui.fbm("CTypes raw"):
            for nr, ct in enumerate(ctypes_generic_list):
                clui.message(str(ct), cont=True)
                sps = lind(spaths, what2what(ctypes_generic, [ct]))
                plot_spaths_traces(fof(sps), name=str(ct) + '_raw', split=False, spp=spp)
                if ct in master_paths:
                    if master_paths[ct] is not None:
                        plot_spaths_traces([master_paths[ct]], name=str(ct) + '_raw_master', split=False, spp=spp,
                                           smooth=lambda anything: anything)

    if options.ctypes_smooth:
        with clui.fbm("CTypes smooth"):
            for nr, ct in enumerate(ctypes_generic_list):
                clui.message(str(ct), cont=True)
                sps = lind(spaths, what2what(ctypes_generic, [ct]))
                plot_spaths_traces(fof(sps), name=str(ct) + '_smooth', split=False, spp=spp, smooth=smooth)
                if ct in master_paths_smooth:
                    if master_paths_smooth[ct] is None: continue
                    plot_spaths_traces([master_paths_smooth[ct]], name=str(ct) + '_smooth_master', split=False, spp=spp,
                                       smooth=lambda anything: anything)
                if ct in master_paths:
                    if master_paths[ct] is None: continue
                    plot_spaths_traces([master_paths[ct]], name=str(ct) + '_raw_master_smooth', split=False, spp=spp,
                                       smooth=smooth)

    fof = lambda sp: list(fractionof(sp, f=make_fracion(options.all_paths_amount,len(sp))))

    if options.all_paths_raw:
        with clui.fbm("All raw paths"):
            plot_spaths_traces(fof(spaths), name='all_raw', split=options.all_paths_split, spp=spp)
    if options.all_paths_raw_io:
        with clui.fbm("All raw paths io"):
            plot_spaths_inlets(fof(spaths), name='all_raw_paths_io', spp=spp)

    if options.all_paths_smooth:
        with clui.fbm("All smooth paths"):
            plot_spaths_traces(fof(spaths), name='all_smooth', split=options.all_paths_split, spp=spp, smooth=smooth)
    if options.all_paths_smooth_io:
        with clui.fbm("All smooth paths io"):
            plot_spaths_inlets(fof(spaths), name='all_smooth_paths_io', spp=spp)

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
