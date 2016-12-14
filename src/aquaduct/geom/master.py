# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2016  Tomasz Magdziarz, Alicja Płuciennik, Michał Stolarczyk <info@aquaduct.pl>
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

# this modlue is a prototype and have to be rewritten

import logging

logger = logging.getLogger(__name__)

import numpy as np
from scipy.spatial.distance import cdist, pdist

from aquaduct.traj.paths import GenericPathTypeCodes, GenericPaths, yield_single_paths, MasterPath
from aquaduct.utils.helpers import list_blocks_to_slices, strech_zip, zip_zip, xzip_xzip
from aquaduct.utils import clui
from aquaduct.utils.maths import make_default_array
from aquaduct.traj.inlets import InletClusterGenericType, InletClusterExtendedType

from multiprocessing import Queue, Manager, Lock, Value, Process

import multiprocessing

from itertools import izip_longest

from functools import partial


def fit_trace_to_points(trace, points):
    dist = cdist(trace, points)
    points_min = np.argmin(dist, 1).tolist()
    points_used = []
    for pm in points_min:
        if pm in points_used:
            continue
        points_used.append(pm)
    return points_used


def decide_on_type(cont, s2o_treshold=0.5):
    # possible types are:
    # GenericPathTypeCodes.object_name
    # GenericPathTypeCodes.scope_name
    o = 0
    if GenericPathTypeCodes.object_name in cont:
        o = cont[GenericPathTypeCodes.object_name]
    s = 0
    if GenericPathTypeCodes.scope_name in cont:
        s = cont[GenericPathTypeCodes.scope_name]
    # get s to o value
    if s == o:
        s2o = 0.5
    else:
        if o == 0:
            s2o = 1.
        else:
            s2o = float(s) / (o + s)
    # decide on type
    # print s,o
    # print s2o, {True:'sco',False:'obj'}[s2o >= s2o_treshold],cont
    if s2o >= s2o_treshold and s > 0:
        return GenericPathTypeCodes.scope_name
    return GenericPathTypeCodes.object_name


def get_weights_(spaths, smooth=None):
    # max len
    max_len = [sp.get_distance_cont(smooth=smooth)[-1] for sp in spaths]
    arg_max_len = np.argmax(max_len)
    max_len = max(max_len)
    # get weights as both lengths
    weights = make_default_array(
        [sz for sz in strech_zip(*[sp.get_distance_both_cont(smooth=smooth, normalize=max_len) for sp in spaths])])
    # add 0.5 - len_both of the lognest path
    length_both_of_max_len = 0.5 - spaths[arg_max_len].get_distance_both_cont(smooth=smooth, normalize=max_len)

    weights = make_default_array([lboml + w for w, lboml in strech_zip(weights, length_both_of_max_len)])
    return weights ** 10


def get_mean_coord_(coords, l):
    # l >> 0
    coord0 = make_default_array(np.median(coords, 0))
    # l >> 1
    coord5 = coords[np.argmax(cdist(coords, np.matrix(coord0)))]

    return make_default_array(np.average([coord0, coord5], 0, [1 - l, l]))


def concatenate(*args):
    for a in args:
        for e in a:
            yield e


part2type_dict = {0: GenericPathTypeCodes.scope_name,
                  1: GenericPathTypeCodes.object_name,
                  2: GenericPathTypeCodes.scope_name}

parts = (0, 1, 2)


class CTypeSpathsCollectionWorker(object):
    # takes group of paths belonging to one ctype and allows to get MasterPath
    def __init__(self, spaths=None, ctype=None, bias_long=5, smooth=None):
        self.spaths = spaths
        assert isinstance(ctype, InletClusterGenericType) or isinstance(ctype, InletClusterExtendedType)
        self.ctype = ctype
        self.bias_long = bias_long
        self.smooth = smooth

        self.lens_cache = None
        self.lens_real_cache = None
        self.lens_norm_cache = None
        self.full_size_cache = None

    def coords_types_prob_widths(self, sp_slices_):
        # get zz coords and types
        coords_zz = [sp.get_coords_cont(smooth=self.smooth)[sl] for sp, sl in zip(self.spaths, sp_slices_)]

        # here we have coords_zz and types_zz
        # and we can calculate coords, types_prob, widths

        # make lens_zz which are lens corrected to the lenghts of coords_zz and normalized to zip_zip number of obejcts
        lens_zz = []
        for l, coord_z in zip(self.lens_cache, coords_zz):
            if len(coord_z) > 0:
                lens_zz.append([float(l) / len(coord_z)] * len(coord_z))  # normalize and correct lengths
            else:
                # lens_zz.append([float(l)] * len(coord_z))
                lens_zz.append([])

        # here we have coords_zz, types_zz, lens_zz
        # and we can calculate coords, types_prob, widths

        # concatenate zip_zip coords and lens
        coords_zz_cat = list(concatenate(*coords_zz))
        del coords_zz

        lens_zz_cat = list(concatenate(*lens_zz))
        del lens_zz
        # average coords_zz_cat using weights of lens_zz_cat

        coords_to_append = make_default_array(np.average(coords_zz_cat, axis=0, weights=lens_zz_cat))
        del lens_zz_cat

        # calculate widths
        if len(self.spaths) > 1:  # is the len of coords_zz the same as sp_slices_ and self.spaths?
            widths_to_append = make_default_array(np.mean(pdist(coords_zz_cat, 'euclidean')))
        else:
            widths_to_append = 0.
        del coords_zz_cat

        # concatenate zip_zip gtypes
        types_zz_cat = list(concatenate(*[sp.gtypes_cont[sl] for sp, sl in zip(self.spaths, sp_slices_)]))
        del sp_slices_
        # append type porbability to types

        types_to_append = float(types_zz_cat.count(GenericPathTypeCodes.scope_name)) / len(types_zz_cat)

        return coords_to_append, types_to_append, widths_to_append

    def __call__(self, nr_sp_slices_):
        return nr_sp_slices_[0], self.coords_types_prob_widths(nr_sp_slices_[-1])


class CTypeSpathsCollection(object):
    parts = (0, 1, 2)  # spath parts

    # takes group of paths belonging to one ctype and allows to get MasterPath
    def __init__(self, spaths=None, ctype=None, bias_long=5, pbar=None, threads=1):

        self.pbar = pbar
        self.threads = threads
        logger.debug("Threads passed %d", threads)

        self.spaths = spaths
        assert isinstance(ctype, InletClusterGenericType) or isinstance(ctype, InletClusterExtendedType)
        self.ctype = ctype
        self.bias_long = bias_long

        # precompute some values
        self.beat()

        with clui.tictoc('spaths props cache in %s' % str(self.ctype)):
            self.lens_cache = self.lens()
            self.lens_real_cache = self.lens_real()
            self.lens_norm_cache = self.lens_norm()
            self.full_size_cache = self.full_size()
        self.beat()

        # self.lock = Lock()

    def beat(self):
        if self.pbar is not None:
            self.pbar.heartbeat()

    def update(self):
        if self.pbar is not None:
            self.pbar.update(1)

    def lens(self):
        # get total lenghts of all paths - needed as weights in averaging
        # if ctype in #:# and not 0 and not None then take object part only length
        if self.ctype.input is not None:
            if self.ctype.input > 0:
                if self.ctype.input == self.ctype.output:
                    return make_default_array([float(len(sp.types_object)) for sp in self.spaths])
        return make_default_array([float(sp.size) for sp in self.spaths])

    def lens_norm(self):
        lens = self.lens()
        if np.max(lens) > 0:
            lens /= np.max(lens)  # normalize
            return lens ** self.bias_long  # bias to long paths

    def lens_real(self):
        return [sp.size for sp in self.spaths]

    def full_size(self):
        # get total size (desired) of master path
        # first check what is the size of paths in all parts and normalize and then scale them
        sizes = []
        for part in self.parts:
            # lengths of all paths of part part
            lens = make_default_array([float(len(sp.types[part])) for sp in self.spaths])
            if np.max(lens) > 0:
                lens /= np.max(lens)  # normalization
                lens = lens ** self.bias_long  # scale them by increasing weights of long paths
            if sum(lens) == 0:
                sizes.append(0)
            else:
                # weighted average by paths lengths
                sizes.append(int(np.average([len(sp.types[part]) for sp in self.spaths], 0, lens)))
        return sum(sizes)  # total size (desired)

    @staticmethod
    def simple_types_distribution(types):
        # possible types are:
        # GenericPathTypeCodes.object_name
        # GenericPathTypeCodes.scope_name
        td_in, td_obj, td_out = 0, 0, 0
        sls = list(list_blocks_to_slices(types))
        if GenericPathTypeCodes.scope_name in types[sls[0]]:
            # this is input part
            td_in = len(types[sls[0]])
        if GenericPathTypeCodes.scope_name in types[sls[-1]]:
            # this is output part
            td_out = len(types[sls[-1]])
        # the rest is object
        td_obj = len(types) - td_in - td_out
        return map(lambda x: float(x) / len(types), (td_in, td_obj, td_out))

    def types_distribution(self):
        # make median distribuitions
        return np.matrix(make_default_array(np.median([self.simple_types_distribution(sp.gtypes_cont) for sp in self.spaths], axis=0)))

    def types_prob_to_types(self, types_prob):
        # get proper types
        types_dist_orig = self.types_distribution()
        types_dist_range = list(set(types_prob))
        types_thresholds = []
        for t in types_dist_range:
            new_pro_types = [{True: GenericPathTypeCodes.scope_name,
                              False: GenericPathTypeCodes.object_name}[typ >= t] for typ in types_prob]
            types_thresholds.append(make_default_array(cdist(np.matrix(self.simple_types_distribution(new_pro_types)),
                                          types_dist_orig, metric='euclidean')))
            self.beat()
        # get threshold for which value of types_thresholds is smallest
        types = [{True: GenericPathTypeCodes.scope_name,
                  False: GenericPathTypeCodes.object_name}[typ >= types_dist_range[np.argmin(types_thresholds)]] for typ
                 in types_prob]
        return types

    def get_master_path(self, smooth=None, resid=0):

        worker = CTypeSpathsCollectionWorker(spaths=self.spaths, ctype=self.ctype, bias_long=self.bias_long,
                                             smooth=smooth)
        # add some spaths properties to worker
        worker.lens_cache = self.lens_cache
        worker.lens_real_cache = self.lens_real_cache
        worker.lens_norm_cache = self.lens_norm_cache
        worker.full_size_cache = self.full_size_cache

        full_size = self.full_size_cache

        # containers for coords, types and widths of master path
        coords = [None] * full_size
        types = [None] * full_size
        widths = [None] * full_size

        # pbar magic
        pbar_previous = 0
        pbar_factor = float(len(self.spaths)) / full_size

        # create pool of workers - mapping function
        map_fun = map
        if self.threads > 1:
            pool = multiprocessing.Pool(self.threads)
            map_fun = pool.imap_unordered
            chunk_size = int(full_size / self.threads ** 2)
            if chunk_size == 0:
                chunk_size = 1
            map_fun = partial(pool.imap_unordered,chunksize=chunk_size)

        # TODO: it is possible to add pbar support here!
        # TODO: consider of using imap_unordered, it might use less memory in this case
        spath_nr_max = 0
        for pbar_nr, (spath_nr, (coords_, types_, widths_)) in enumerate(
                map_fun(worker, enumerate(xzip_xzip(*worker.lens_real_cache, N=full_size)))):
            coords[spath_nr] = coords_
            types[spath_nr] = types_
            widths[spath_nr] = widths_
            spath_nr_max = max(spath_nr, spath_nr_max)
            pbar_current = int((pbar_nr + 1) * pbar_factor)
            if pbar_current > pbar_previous:
                pbar_previous = pbar_current
                self.update()  # update progress bar
            else:
                self.beat()
        assert pbar_nr == spath_nr_max, "Internal error. Final global progress of master path generation not synced with maximal number of spath. Please send a bug report to developer(s): %s" % clui.mail

        if self.threads > 1:
            pool.close()
            pool.join()
            pool.terminate()
            del pool

        # at this stage we have coords, widths and types probability

        # get proper types
        with clui.tictoc('proper tests in %s' % str(self.ctype)):
            types = self.types_prob_to_types(types)

        # make frames
        frames = range(len(coords))

        # finalize

        # max min frames
        min_pf = 0
        max_pf = len(coords) - 1
        if self.ctype is None:  # this never happens because of assertion in __init__
            min_pf = None
            max_pf = None
        else:
            if self.ctype.input is not None:
                min_pf = None
            if self.ctype.output is not None:
                max_pf = None

        with clui.tictoc('generic paths in %s' % str(self.ctype)):
            # get and populate GenericPath
            gp = GenericPaths(resid, min_pf=min_pf, max_pf=max_pf)
            for c, t, f in zip(coords, types, frames):  # TODO: remove loop
                gp.add_type(f, t)
                gp.add_coord(c)
        # now try to get first SinglePath, if unable issue WARNING
        with clui.tictoc('separate paths in %s' % str(self.ctype)):
            try:
                sp = list(yield_single_paths([gp]))[0]
            except IndexError:
                logger.warning('No master path found for ctype %s' % str(self.ctype))
                return None
        # finally get MasterPath and add widths
        mp = MasterPath(sp)
        mp.add_width(widths)
        return mp


def create_master_spath(spaths, smooth=None, resid=0, ctype=None, bias_long=5, pbar=None):
    def beat():
        if pbar is not None:
            pbar.heartbeat()

    def update():
        if pbar is not None:
            pbar.update(1)

    # first check what is the size of paths in all parts and normalize and then scale them
    sizes = []
    for part in parts:
        # lengths of all paths of part part
        lens = make_default_array([float(len(sp.types[part])) for sp in spaths])
        if np.max(lens) > 0:
            lens /= np.max(lens)  # normalization
            lens = lens ** bias_long  # scale them by increasing weights of long paths
        if sum(lens) == 0:
            sizes.append(0)
        else:
            # weighted average by paths lengths
            sizes.append(int(np.average([len(sp.types[part]) for sp in spaths], 0, lens)))
        beat()
    full_size = sum(sizes)  # total size (desired)
    pbar_factor = float(len(spaths)) / full_size

    # get total lenghts of all paths - needed as weights in averaging
    lens = make_default_array([float(sp.size) for sp in spaths])
    # if ctype in #:# and not 0 and not None then take object part only length
    if ctype is not None:
        if ctype.input is not None:
            if ctype.input > 0:
                if ctype.input == ctype.output:
                    lens = make_default_array([float(len(sp.types_object)) for sp in spaths])
    # normalize and incearse weight of long paths (depends ob ctype)
    if np.max(lens) > 0:
        lens /= np.max(lens)  # normalize
        lens = lens ** bias_long  # bias to long paths

    # containers for coords, types and widths of master path
    coords = []
    types = []
    widths = []
    pbar_previous = 0

    # following loop calculates types distribution only
    # loop over zip zipped [smooth] coords of all paths and gtypes with size set to full_size
    for pbar_nr, (coords_zz, types_zz) in enumerate(
            zip(zip_zip(*[sp.get_coords_cont(smooth=smooth) for sp in spaths], N=full_size),
                zip_zip(*[sp.gtypes_cont for sp in spaths], N=full_size))):
        # make lens_zz which are lens corrected to the lenghts of coords_zz and normalized to zip_zip number of obejcts
        lens_zz = []
        for l, coord_z in zip(lens, coords_zz):
            if len(coord_z) > 0:
                lens_zz.append([float(l) / len(coord_z)] * len(coord_z))  # normalize and correct lengths
            else:
                # lens_zz.append([float(l)] * len(coord_z))
                lens_zz.append([])
        # concatenate zip_zip coords and lens
        coords_zz_cat = list(concatenate(*coords_zz))
        lens_zz_cat = list(concatenate(*lens_zz))
        # average coords_zz_cat using weights of lens_zz_cat
        coords.append(make_default_array(np.average(coords_zz_cat, 0, lens_zz_cat)))
        # calculate widths
        if len(coords_zz) > 1:
            # try tu use weighted distance - wminkowski with p=2 is equivalent to weighted euclidean
            # id_of_max = np.argmax(pdist(coords_zz_cat, 'wminkowski', p=2, w=lens_zz_cat))
            # widths.append(pdist(coords_zz_cat, 'euclidean')[id_of_max])
            widths.append(np.mean(pdist(coords_zz_cat, 'euclidean')))
        else:
            widths.append(0.)
        # concatenate zip_zip gtypes
        types_zz_cat = list(concatenate(*types_zz))
        # # pick correct type..., check distance of coords[-1] to coords_zz_cat
        # types_cdist = cdist(np.matrix(coords[-1]), coords_zz_cat, metric='euclidean')
        # types.append(types_zz_cat[np.argmin(types_cdist)])
        # types.append(decide_on_type(Counter(types_zz_cat)))
        types.append(float(types_zz_cat.count(GenericPathTypeCodes.scope_name)) / len(types_zz_cat))

        pbar_current = int((pbar_nr + 1) * pbar_factor)
        if pbar_current > pbar_previous:
            pbar_previous = pbar_current
            update()  # update progress bar
        else:
            beat()
    # get proper types
    # make median distribuitions
    types_dist_orig = np.matrix(
        make_default_array(np.median([CTypeSpathsCollection.simple_types_distribution(sp.gtypes_cont) for sp in spaths], axis=0)))
    types_dist_range = list(set(types))
    types_thresholds = []
    for t in types_dist_range:
        new_pro_types = [{True: GenericPathTypeCodes.scope_name,
                          False: GenericPathTypeCodes.object_name}[typ >= t] for typ in types]
        types_thresholds.append(make_default_array(cdist(np.matrix(CTypeSpathsCollection.simple_types_distribution(new_pro_types)),
                                      types_dist_orig, metric='euclidean')))
        beat()
    # get threshold for which value of types_thresholds is smallest
    types = [{True: GenericPathTypeCodes.scope_name,
              False: GenericPathTypeCodes.object_name}[typ >= types_dist_range[np.argmin(types_thresholds)]] for typ in
             types]

    frames = range(len(coords))

    min_pf = 0
    max_pf = len(coords) - 1

    if ctype is None:
        min_pf = None
        max_pf = None
    else:
        if ctype.input is not None:
            min_pf = None
        if ctype.output is not None:
            max_pf = None

    gp = GenericPaths(resid, min_pf=min_pf, max_pf=max_pf)
    for c, t, f in zip(coords, types, frames):  # TODO: remove loop
        gp.add_type(f, t)
        gp.add_coord(c)
        beat()
    # return gp
    beat()  # touch progress bar
    try:
        sp = list(yield_single_paths([gp]))[0]
    except IndexError:
        logger.warning('No master path found for ctype %s' % str(ctype))
        return None
    beat()  # touch progress bar
    mp = MasterPath(sp)
    mp.add_width(widths)
    return mp


def calculate_master(spaths_resid_ctype_smooth):
    spaths, resid, ctype, smooth = spaths_resid_ctype_smooth
    return CTypeSpathsCollection(spaths=spaths, ctype=ctype).get_master_path(smooth=smooth, resid=resid)
