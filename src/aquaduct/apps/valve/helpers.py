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

import operator
import os
import re
import sys
import numpy as np

from aquaduct.apps.data import CRIC
from aquaduct import logger
from aquaduct.apps.valve import asep
from aquaduct.geom import traces
from aquaduct.geom.cluster_available_methods import AVAILABLE_METHODS as available_clustering_methods
from aquaduct.geom.cluster import BarberCluster, PerformClustering
from aquaduct.geom.smooth import WindowSmooth, ActiveWindowSmooth, DistanceWindowSmooth, MaxStepSmooth, \
    WindowOverMaxStepSmooth, ActiveWindowOverMaxStepSmooth, DistanceWindowOverMaxStepSmooth, SavgolSmooth
from sklearn.cluster import DBSCAN, AffinityPropagation, KMeans, MeanShift, Birch
from aquaduct.traj.barber import WhereToCut
from aquaduct.utils import clui
from aquaduct.utils.helpers import is_number, Auto


class NP(object):
    def __init__(self, pbar, next_len=False):
        self.paths = list()
        self.pbar = pbar
        self.next_len = next_len

    def reinit(self, pbar, next_len=False):
        self.pbar = pbar
        self.next_len = next_len

    def next(self, n):
        if n > 1 and self.next_len:
            self.pbar.next(step=n)
        else:
            self.pbar.next()

    def callback_cric_next(self, result):
        CRIC.update_cric(result[-1])
        self.paths.extend(result[:-1])
        self.next(len(result) - 1)

    def callback_next(self, result):
        self.paths.extend(result)
        self.next(len(result))

    def callback_append_next(self, result):
        self.paths.append(result)
        self.pbar.next()


def get_res_in_scope(is_res_in_scope, res):
    res_new = None
    for iris, r in zip(is_res_in_scope, res.iterate_over_residues()):
        if iris:
            if res_new is None:
                res_new = r
            else:
                res_new += r
    return res_new


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
    assert coptions.method in available_clustering_methods, 'Unknown clustering method %s.' % coptions.method

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
        assert isinstance(loption, str), "Wrong Linearize method definition: %r" % loption
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


def discard_short_etc(spaths, short_paths=None, short_object=None, short_logic=None):
    # worker for discarding paths, can be time consuming because of filling coords cache
    if short_object is not None:
        return [sp for sp in spaths if
                short_logic(sp.size > short_paths, sp.object_len > short_object)] + [CRIC]
        # return len(spaths), [sp for sp in spaths if
        #                     short_logic(sp.size > short_paths, sp.object_len > short_object)], CRIC
    else:
        return [sp for sp in spaths if sp.size > short_paths]
        # return len(spaths), [sp for sp in spaths if sp.size > short_paths]


def center_of_object(spath):
    # import ipdb as pdb; pdb.set_trace()
    return 1, spath.center_of_object, CRIC


def get_allow_size_function(rt=None):
    if not isinstance(rt, str):
        return None
    assert re.compile('^[<>=]+[0-9.]+$').match(rt) is not None, "Wrong threshold definition: %s" % rt
    op = re.compile('[<>=]+')
    op = ''.join(sorted(op.findall(rt)[0]))
    vl = re.compile('[0-9.]+')
    vl = float(vl.findall(rt)[0])
    operator_dict = {'>': operator.gt,
                     '=>': operator.ge,
                     '<=': operator.le,
                     '<': operator.lt}
    assert op in list(operator_dict.keys()), "Unsupported operator %s in threshold %s" % (op, rt)
    return lambda size_of_cluster: operator_dict[op](size_of_cluster, vl)


def potentially_recursive_clustering(config=None,
                                     clustering_name=None,
                                     inlets_object=None,
                                     spaths=None,
                                     message='clustering',
                                     deep=0,
                                     max_level=5):
    with clui.fbm("Performing %s, level %d of %d" % (message, deep, max_level), cont=False):
        logger.debug('Clustering options section: %s' % clustering_name)
        cluster_options = config.get_cluster_options(section_name=clustering_name)
        clui.message('Clustering options:')
        for k, v in cluster_options._asdict().items():
            clui.message("%s = %s" % (str(k), str(v)))
        # TODO: Print clustering options in a nice way!
        clustering_function = get_clustering_method(cluster_options, config)
        # special case of barber!!!
        if cluster_options.method == 'barber':
            logger.debug('Getting inltets refs...')
            inlets_refs = set(inlets_object.get_inlets_references())
            logger.debug('Starting wtc...')
            wtc = WhereToCut(spaths=(sp for sp in spaths if sp.id in inlets_refs),
                             expected_nr_of_spaths=len(inlets_refs),
                             forceempty=True,
                             **clustering_function.method_kwargs)
            # clouds = wtc.cloud_groups(progress=True)
            logger.debug('Getting spheres...')
            inlets_object.add_spheres(wtc.spheres)
        logger.debug('Proceed with clustering, skip size...')
        # get skip_size function according to recursive_treshold
        skip_size = SkipSizeFunction(cluster_options.recursive_threshold)
        logger.debug('Proceed with clustering, call method...')
        inlets_object.perform_reclustering(clustering_function, skip_outliers=True, skip_size=skip_size)
    clui.message('Number of clusters detected so far: %d' % len(inlets_object.clusters_list))

    if cluster_options.recursive_clustering:
        deep += 1
        if deep > max_level:
            return
        return potentially_recursive_clustering(config=config,
                                                clustering_name=cluster_options.recursive_clustering,
                                                inlets_object=inlets_object,
                                                spaths=spaths,
                                                deep=deep,
                                                max_level=max_level)


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


class SkipSizeFunction(object):

    def __init__(self, ths_def):

        self.thresholds = []
        if isinstance(ths_def, str):
            for thd in ths_def.split():
                self.thresholds.append(get_allow_size_function(thd))

    def __call__(self, size_of_cluster):
        for thd in self.thresholds:
            if not thd(size_of_cluster):
                return True
        return False


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
            info2print = (self.nr_template % (nr + 1)) + info2print
        if self.output2stderr:
            clui.message(info2print)
        print(info2print, file=self.filehandle)

    def sep(self):
        self(asep())

    def thead(self, info2print):
        self(clui.thead(info2print))

    def tend(self, info2print):
        self(clui.tsep(info2print))

    def under(self, info2print):
        self(clui.underline(info2print))


def results_n(rn):
    if isinstance(rn, np.ndarray):
        return rn
    else:
        return np.memmap(rn[0], mode='r', dtype=np.int8, shape=rn[1])
