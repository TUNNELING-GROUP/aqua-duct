# -*- coding: utf-8 -*-

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

import gc
from itertools import chain

import numpy as np
from aquaduct.apps.data import GCS
from aquaduct import logger

from setuptools import find_packages, setup
from aquaduct.traj.sandwich import open_traj_reader, ResidueSelection
from aquaduct.utils.helpers import create_tmpfile, iterate_or_die

#from aquaduct.apps.valve.core import GenericPaths
from aquaduct.traj.paths import GenericPaths
from aquaduct.apps.valve.helpers import results_n


class assign_nonsandwiched_paths(object):
    """
    Worker which assign non-sandwiched paths to object container
    """

    def __call__(self, args):
        """
        :param args: tuple with object and paths
        :return: aquaduct.traj.paths.GenericPaths
        """
        pat, nfos = args
        pat.add_012(nfos)
        return pat


class assign_sandwiched_paths(object):
    """
    Worker which assign sandwiched paths to object container
    """

    def __init__(self, all_res_ids, all_res_names, max_pf, results):
        """
        Constructor
        :param all_res_ids: residues ids
        :param all_res_names: residues names
        :param max_pf: maximum possible frame
        :param results: residue coords
        :param pbar: progress bar
        """
        self.all_res_ids = all_res_ids
        self.all_res_names = all_res_names
        self.max_pf = max_pf
        self.results = results

    def __call__(self, pnr):
        """
        :param pnr: residue id
        :return: aquaduct.traj.paths.GenericPaths
        """
        new_p = GenericPaths((0, self.all_res_ids[pnr]),
                             name_of_res=self.all_res_names[pnr],
                             min_pf=0, max_pf=self.max_pf)

        new_p.add_012(
            np.fromiter(
                chain(*(results_n(self.results[n])[:, pnr] for n in sorted(self.results.keys()))),
                dtype=np.int8))

        return new_p


def stage_I_worker_q(input_queue, results_queue, pbar_queue):
    for input_data in iter(input_queue.get, None):
        layer_number, traj_reader_proto, scope_everyframe, scope, scope_convexhull, scope_convexhull_inflate, object_selection, progress_freq = input_data

        center_of_system = np.zeros(3)
        all_res = None

        traj_reader = open_traj_reader(traj_reader_proto)

        if not scope_everyframe:
            scope_selection = traj_reader.parse_selection(scope)

        progress = 0
        progress_freq_flex = min(1, progress_freq)
        frame_rid_in_object = []
        # the loop over frames
        for frame in traj_reader.iterate_over_frames():
            if scope_everyframe:
                scope_selection = traj_reader.parse_selection(scope)
            # center of system
            center_of_system += scope_selection.center_of_mass()
            # current res selection
            res = traj_reader.parse_selection(object_selection).residues()
            # find matching residues, ie those which are in the scope:
            res_new = scope_selection.containing_residues(res, convex_hull=scope_convexhull,
                                                convex_hull_inflate=scope_convexhull_inflate)
            res_new.uniquify()  # here is a list of residues in this layer that are in the object and in the scope
            # TODO: change way of center_of_system calculation to reflect center of object?
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
            progress += 1
            if progress == progress_freq_flex:
                pbar_queue.put(progress)
                progress = 0
                progress_freq_flex = min(progress_freq_flex * 2, progress_freq)
        # sent results
        results_queue.put({layer_number: (all_res, frame_rid_in_object, center_of_system)})
        pbar_queue.put(progress)


def stage_II_worker_q(input_queue, results_queue, pbar_queue):
    # input queue loop
    for input_data in iter(input_queue.get, None):
        # get data
        layer_number, traj_reader_proto, scope_everyframe, scope, scope_convexhull, scope_convexhull_inflate, object_selection, all_res_this_layer, frame_rid_in_object, is_number_frame_rid_in_object, progress_freq = input_data
        # open trajectory and get number of frames
        traj_reader = open_traj_reader(traj_reader_proto)
        number_of_frames = traj_reader.number_of_frames()
        # scope is evaluated only once before loop over frames so it cannot be frame dependent
        if not scope_everyframe:
            scope = traj_reader.parse_selection(scope)
            logger.debug("Scope definition evaluated only once for given layer")
        else:
            logger.debug("Scope definition evaluated in every frame, this might be very slow.")
        # get ids only once
        all_res_this_ids = list(all_res_this_layer.ids())
        # big container for 012 path data
        # cache???
        if GCS.cachedir:
            number_frame_object_scope = np.memmap(create_tmpfile(ext='dat', dir=GCS.cachedir),
                                                  dtype=np.int8,
                                                  shape=(number_of_frames, all_res_this_layer.len()))
        else:
            number_frame_object_scope = np.zeros((number_of_frames, all_res_this_layer.len()),
                                                 dtype=np.int8)
        # progress reported peridicaly
        progress = 0
        progress_gc = 0
        progress_freq_flex = min(1, progress_freq)
        # the loop over frames, use izip otherwise iteration over frames does not work
        for rid_in_object, frame in zip(iterate_or_die(frame_rid_in_object, times=number_of_frames),
                                         traj_reader.iterate_over_frames()):
            # do we have object data?
            if not is_number_frame_rid_in_object:
                rid_in_object = [rid[-1] for rid in traj_reader.parse_selection(object_selection).residues().ids()]
            # assert rid_in_object is not None
            is_res_in_object = (rid[-1] in rid_in_object for rid in all_res_this_ids)
            # should scope be evaluated?
            if scope_everyframe:
                scope = traj_reader.parse_selection(scope)
            # check if all_res are in the scope, reuse res_ids_in_object_over_frames
            is_res_in_scope = scope.contains_residues(all_res_this_layer, convex_hull=scope_convexhull,
                                                      convex_hull_inflate=scope_convexhull_inflate,
                                                      known_true=None)  # known_true could be rid_in_object
            # store results in the container
            number_frame_object_scope[frame, :] = np.array(list(map(sum, zip(is_res_in_object, is_res_in_scope))),
                                                           dtype=np.int8)
            # increase progress counter and report progress if needed
            progress += 1
            progress_gc += 1
            if progress == progress_freq_flex:
                pbar_queue.put(progress)
                progress = 0
                progress_freq_flex = min(progress_freq_flex * 2, progress_freq)
                if GCS.cachedir:
                    number_frame_object_scope.flush()
            if progress_gc == progress_freq * 100:
                gc.collect()
                progress_gc = 0
        # sent results to results_queue
        # cache???
        if GCS.cachedir:
            number_frame_object_scope.flush()
            results_queue.put({layer_number: (number_frame_object_scope.filename, number_frame_object_scope.shape)})
            del number_frame_object_scope
        else:
            results_queue.put({layer_number: number_frame_object_scope})
        if progress:
            pbar_queue.put(progress)
        # termination
        gc.collect()


def stage_II_worker_q_twoways(input_queue, results_queue, pbar_queue):
    # input queue loop
    for input_data in iter(input_queue.get, None):
        # get data
        layer_number, traj_reader_proto, scope_everyframe, scope, scope_convexhull, scope_convexhull_inflate, object_selection, all_res_this_layer, frame_rid_in_object, is_number_frame_rid_in_object, progress_freq = input_data
        # open trajectory and get number of frames
        traj_reader = open_traj_reader(traj_reader_proto)
        number_of_frames = traj_reader.number_of_frames()
        # scope is evaluated only once before loop over frames so it cannot be frame dependent
        if not scope_everyframe:
            scope = traj_reader.parse_selection(scope)
            logger.debug("Scope definition evaluated only once for given layer")
        else:
            logger.debug("Scope definition evaluated in every frame, this might be very slow.")
        # get ids only once
        all_res_this_ids = list(all_res_this_layer.ids())
        # big container for 012 path data
        # cache???
        if GCS.cachedir:
            number_frame_object_scope = np.memmap(create_tmpfile(ext='dat', dir=GCS.cachedir),
                                                  dtype=np.int8,
                                                  shape=(number_of_frames, all_res_this_layer.len()))
        else:
            number_frame_object_scope = np.zeros((number_of_frames, all_res_this_layer.len()),
                                                 dtype=np.int8)
        # the loop over frames, use izip otherwise iteration over frames does not work
        progress = 0
        for reverse in (False, True):
            # progress reported peridicaly
            progress_gc = 0
            progress_freq_flex = min(1, progress_freq)
            # which all_res should be evalated?
            all_res_eval = np.ones(len(all_res_this_ids), dtype=bool)
            for rid_in_object, frame in zip(
                    iterate_or_die(frame_rid_in_object, times=number_of_frames, reverse=reverse),
                    traj_reader.iterate_over_frames(reverse=reverse)):
                # do we have object data?
                if not is_number_frame_rid_in_object:
                    rid_in_object = [rid[-1] for rid in traj_reader.parse_selection(object_selection).residues().ids()]
                # assert rid_in_object is not None
                is_res_in_object = np.fromiter((rid[-1] in rid_in_object for rid in all_res_this_ids), dtype=bool)
                # if in object do not do scope check
                all_res_eval[is_res_in_object] = False
                if reverse:
                    all_res_eval[number_frame_object_scope[frame, :] > 0] = False
                all_res_this_ids_eval = (i[-1] for te, i in zip(all_res_eval, all_res_this_ids) if te)
                all_res_this_layer_eval = ResidueSelection(
                    {all_res_this_layer.numbers()[0]: list(all_res_this_ids_eval)})
                # should scope be evaluated?
                if scope_everyframe:
                    scope = traj_reader.parse_selection(scope)
                # check if all_res are in the scope, reuse res_ids_in_object_over_frames
                is_res_in_scope_eval = scope.contains_residues(all_res_this_layer_eval, convex_hull=scope_convexhull,
                                                               convex_hull_inflate=scope_convexhull_inflate,
                                                               known_true=None)  # known_true could be rid_in_object
                is_res_in_scope = np.zeros(len(all_res_this_ids), dtype=bool)
                is_res_in_scope[
                    [int(nr) for nr, iris in zip(np.argwhere(all_res_eval), is_res_in_scope_eval) if iris]] = True
                if not reverse:
                    is_res_in_scope[is_res_in_object] = True
                # store results in the container
                if reverse:
                    number_frame_object_scope[frame, is_res_in_scope] = 1
                    is_res_in_scope[is_res_in_object] = True
                else:
                    number_frame_object_scope[frame, :] = np.array(list(map(sum, zip(is_res_in_object, is_res_in_scope))),
                                                                   dtype=np.int8)
                # increase progress counter and report progress if needed
                progress += 1
                progress_gc += 1
                if progress == progress_freq_flex:
                    pbar_queue.put(progress * 0.5)
                    progress = 0
                    progress_freq_flex = min(progress_freq_flex * 2, progress_freq)
                    if GCS.cachedir:
                        number_frame_object_scope.flush()
                if progress_gc > progress_freq * 100:
                    gc.collect()
                    progress_gc = 0
                all_res_eval = np.zeros(len(all_res_this_ids), dtype=bool)
                all_res_eval[is_res_in_scope] = True

        # sent results to results_queue
        # cache???
        if GCS.cachedir:
            number_frame_object_scope.flush()
            results_queue.put({layer_number: (number_frame_object_scope.filename, number_frame_object_scope.shape)})
            del number_frame_object_scope
        else:
            results_queue.put({layer_number: number_frame_object_scope})
        gc.collect()
        pbar_queue.put(progress * 0.5)
        # termination
