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

import gc
from functools import partial
from itertools import zip_longest, chain
from multiprocessing import Pool, Queue, Process

from scipy.spatial.distance import cdist

from aquaduct import version_nice as aquaduct_version_nice
from aquaduct.apps.data import GCS, save_cric
from aquaduct.apps.valve import sep
from aquaduct.apps.valve.data import get_vda_reader
from aquaduct.geom.master import CTypeSpathsCollection
from aquaduct.traj.barber import barber_paths
from aquaduct.traj.inlets import InletClusterGenericType
from aquaduct.traj.inlets import Inlets
from aquaduct.traj.paths import GenericPaths, yield_single_paths, SinglePath, MacroMolPath
from aquaduct.traj.paths import yield_generic_paths, correct_spaths_ids
from aquaduct.traj.sandwich import ResidueSelection, Reader, mda_ver




from aquaduct.utils.clui import roman
from aquaduct.utils.helpers import iterate_or_die, fractionof, make_fractionof, make_fraction
from aquaduct.utils.helpers import range2int, what2what, lind, robust_and, robust_or
from aquaduct.utils.multip import optimal_threads
from aquaduct.utils.sets import intersection_full

from aquaduct.apps.valve.worker import stage_I_worker_q, stage_II_worker_q, stage_II_worker_q_twoways, \
    assign_nonsandwiched_paths, assign_sandwiched_paths
from aquaduct.apps.valve.helpers import *
from aquaduct.apps.valve.spath import *
from aquaduct.apps.valve.clusters import *
from aquaduct.visual import cmaps

__mail__ = 'info@aquaduct.pl'
__version__ = aquaduct_version_nice()


################################################################################

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
                return dict(((key, val) for key, val in result.items() if 'options' not in key))


################################################################################
# asorted helpers

def get_traced_names(some_paths):
    # some_paths might me paths of spaths or both
    return tuple(sorted(list(set((sp.id.name if isinstance(sp, MacroMolPath) else sp.name for sp in some_paths)))))


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
            range(optimal_threads.threads_count)]
    [p.start() for p in pool]

    if options.scope_convexhull_inflate:
        logger.info("Inflate convex hull by %0.1f" % float(options.scope_convexhull_inflate))
    # feed input_queue with data
    for results_count, (number, traj_reader) in enumerate(Reader.iterate(number=True)):
        input_queue.put((number, traj_reader,
                         options.scope_everyframe, options.scope, options.scope_convexhull,
                         options.scope_convexhull_inflate,
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

    unsandwichize = not Reader.sandwich_mode and len(number_frame_rid_in_object) > 1
    if unsandwichize:
        with clui.fbm("Unsandwichize traced residues"):
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

def waterfall_me(paths, pbar=None):


    number_of_frames = Reader.number_of_frames(onelayer=True)


    N = len(paths)  # number of paths
    for n in range(N):
        frames = paths[0].frames
        types = paths[0].types
        min_pf = 0
        # print "initial",paths[0].id, (paths[0].min_possible_frame, paths[0].max_possible_frame),(frames[0],frames[-1])
        if len(frames):
            for max_pf in list(Reader.edges) + [number_of_frames - 1]:
                # min_pf and max_pf are limits of frames
                this_frames = []
                this_types = []
                for f, t in zip(frames, types):
                    if f >= min_pf and f <= max_pf:
                        this_frames.append(f)
                        this_types.append(t)
                if len(this_frames):
                    #p = GenericPaths(paths[0].id,
                    #                 name_of_res=paths[0].name,
                    #                 min_pf=0,
                    #                 max_pf=number_of_frames-1)
                    p = GenericPaths(paths[0].id,
                                     name_of_res=paths[0].name,
                                     min_pf=min_pf,
                                     max_pf=max_pf)
                    p.update_types_frames(this_types, this_frames)
                    # print p.min_possible_frame,p.max_possible_frame,p.frames[0],p.frames[-1]
                    paths.append(p)
                min_pf = max_pf + 1
        if pbar:
            pbar.next()
        paths.pop(0)


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
            allow_twoway = (not Reader.sandwich_mode) and (sI.scope == sII.scope) and (
                    sI.scope_convexhull == sII.scope_convexhull) and (sI.object == sII.object) and (
                                   sIII.allow_passing_paths == False)

        # prepare and start pool of workers
        if allow_twoway:
            logger.info('Twoway trajectory scan enabled.')
            pool = [Process(target=stage_II_worker_q_twoways, args=(input_queue, results_queue, pbar_queue)) for dummy
                    in
                    range(optimal_threads.threads_count)]
        else:
            pool = [Process(target=stage_II_worker_q, args=(input_queue, results_queue, pbar_queue)) for dummy in
                    range(optimal_threads.threads_count)]
        [p.start() for p in pool]

        if options.scope_convexhull_inflate:
            logger.info("Inflate convex hull by %0.1f" % float(options.scope_convexhull_inflate))
        # feed input_queue with data
        if Reader.sandwich_mode:
            for results_count, (frame_rid_in_object, (number, traj_reader)) in enumerate(
                    zip(iterate_or_die(number_frame_rid_in_object,
                                        times=Reader.number_of_layers()),
                         Reader.iterate(number=True))):
                all_res_layer = all_res.layer(number)

                input_queue.put((number, traj_reader,
                                 options.scope_everyframe, options.scope, options.scope_convexhull,
                                 options.scope_convexhull_inflate,
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
                                 options.scope_everyframe, options.scope, options.scope_convexhull,
                                 options.scope_convexhull_inflate,
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
            if progress % (max(1, optimal_threads.threads_count) ** 2 * 1000) == 0:
                gc.collect()

        # [stop workers]
        [input_queue.put(None) for p in pool]

    # collect results
    pbar = clui.pbar((results_count + 1) + 1, 'Collecting results from layers:')
    results = {}
    for nr, result in enumerate(iter(results_queue.get, None)):
        for rk in list(result.keys()):
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

    unsandwichize = not Reader.sandwich_mode and len(results) > 1
    if not unsandwichize:
        with clui.pbar(len(all_res), 'Creating raw paths:') as pbar:
            new_paths = NP(pbar)
            for number in sorted(results.keys()):
                all_res_layer = all_res.layer(number)

                paths_this_layer = (GenericPaths(resid,
                                                 name_of_res=resname,
                                                 min_pf=0, max_pf=number_of_frames - 1)
                                    for resid, resname in zip(all_res_layer.ids(),
                                                               all_res_layer.names()))

                pool = Pool(processes=optimal_threads.threads_count)
                #list(map(new_paths.callback_append_next, map(assign_nonsandwiched_paths(),
                #                                                             zip(paths_this_layer,
                #                                                                 results_n(results[number]).T))))
                list(map(new_paths.callback_append_next, pool.imap_unordered(assign_nonsandwiched_paths(),
                                                                             zip(paths_this_layer,
                                                                                 results_n(results[number]).T))))

                pool.close()
                pool.join()

    elif unsandwichize:
        # make coherent paths
        # paths names, paths
        max_pf = Reader.number_of_frames() - 1
        # frames_offset = np.cumsum([0] + frames_offset).tolist()[:len(numbers)]

        # pbar = clui.pbar(len(results[numbers[0]]), 'Sandwich deconvolution:')
        with clui.pbar(len(all_res_ids), 'Creating raw paths (sandwich deconvolution):') as pbar:
            new_paths = NP(pbar)
            pool_func = assign_sandwiched_paths(all_res_ids, all_res_names, max_pf, results)

            pool = Pool(processes=optimal_threads.threads_count)
            list(map(new_paths.callback_append_next, pool.imap_unordered(pool_func, list(range(len(all_res_ids))))))

            pool.close()
            pool.join()
    paths = new_paths.paths

    # if edges then split each of path into smaller parts and adjust min_pf and max_pf
    if Reader.edges:
        with clui.pbar(len(paths), 'Waterfall fall:') as pbar:
            waterfall_me(paths, pbar)

    # rm tmp files
    for rn in results.values():
        if not isinstance(rn, np.ndarray):
            os.unlink(rn[0])

    if options.discard_singletons:
        with clui.pbar(maxval=len(paths),
                       mess="Discard singletons (%d) paths:" % int(options.discard_singletons)) as pbar:
            for pat in paths:
                pat.discard_singletons(singl=int(options.discard_singletons))
                pbar.next()

    if options.discard_empty_paths:
        with clui.pbar(maxval=len(paths), mess="Discard residues with empty paths:") as pbar:
            paths = [pat for pat in paths if len(pat.frames) > 0 if pbar.next() is None]

    clui.message("Number of paths: %d" % len(paths))

    return {'all_res': all_res, 'paths': paths, 'options': options._asdict()}


################################################################################

def spaths_paths_ids(spaths):
    seen = []
    for sp in spaths:
        if sp.id.id not in seen:
            seen.append(sp.id.id)
            yield sp.id.id


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
        for sps_nrs in map(ysp, (paths[i:i + n] for i in range(0, len(paths), n))):
            # for sps_nrs in pool.imap_unordered(ysp, (paths[i:i + n] for i in xrange(0, len(paths), n))):
            nr = 0  # if no spaths were returned
            for sp, nr in sps_nrs:
                spaths.append(sp)
                pbar.update(nr_all + nr + 1)
            nr_all += nr + 1
        pool.close()
        pool.join()
        gc.collect()

    if Reader.edges:
        with clui.pbar(len(spaths), "Clean IDs:") as pbar:
            correct_spaths_ids(spaths, pbar)

    clui.message("Created %d separate paths out of %d raw paths" %
                 (len(spaths), len(paths)))

    ######################################################################

    if options.remove_unused_parts:
        with clui.pbar(len(spaths) + len(paths), "Removing unused parts of paths:") as pbar:
            paths = yield_generic_paths(spaths, progress=pbar)
        if Reader.edges:
            with clui.pbar(len(paths), 'Waterfall fall:') as pbar:
                waterfall_me(paths, pbar)

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
                new_spaths = NP(None)
                new_spaths.reinit(pbar)
                # CRIC AWARE MP!
                if short_object is not None:
                    Reader.reset()
                    for pid in list(spaths_paths_ids(spaths)):
                        pool.apply_async(dse, args=(
                            [spaths.pop(nr) for nr in range(len(spaths) - 1, -1, -1) if spaths[nr].id.id == pid],),
                                         callback=new_spaths.callback_cric_next)
                else:
                    for pid in list(spaths_paths_ids(spaths)):
                        pool.apply_async(dse, args=(
                            [spaths.pop(nr) for nr in range(len(spaths) - 1, -1, -1) if spaths[nr].id.id == pid],),
                                         callback=new_spaths.callback_next)
                pool.close()
                pool.join()
                save_cric()
                gc.collect()
                spaths = new_spaths.paths
                del new_spaths
            spaths_nr_new = len(spaths)
            if spaths_nr == spaths_nr_new:
                clui.message("No paths were discarded.")
            else:
                clui.message("%d paths were discarded." % (spaths_nr - spaths_nr_new))
                if options.remove_unused_parts:
                    with clui.pbar(len(spaths) + len(paths), "Removing (again) unused parts of paths:") as pbar:
                        paths = yield_generic_paths(spaths, progress=pbar)
                    if Reader.edges:
                        with clui.pbar(len(paths), 'Waterfall fall:') as pbar:
                            waterfall_me(paths, pbar)
        else:
            clui.message("No paths were discarded - no values were set.")

    ######################################################################

    traced_names = get_traced_names(spaths)

    def iter_over_tn():
        if options.separate_barber:
            for tn in traced_names:
                yield [tn]
        else:
            yield traced_names

    if options.auto_barber:
        new_paths = NP(None)  # no progress bar (yet)
        for tn in iter_over_tn():  # tn is a list of traced names
            with clui.fbm("AutoBarber calculations for %s" % ' '.join(tn), cont=False):
                wtc = WhereToCut(spaths=[sp for sp in spaths if sp.id.name in tn],
                                 selection=options.auto_barber,
                                 mincut=options.auto_barber_mincut,
                                 mincut_level=options.auto_barber_mincut_level,
                                 maxcut=options.auto_barber_maxcut,
                                 maxcut_level=options.auto_barber_maxcut_level,
                                 tovdw=options.auto_barber_tovdw)
                # cut thyself!
                wtc.cut_thyself()  # wtc for given tn

                # how many paths of tn we have?
                tnpaths = [nr for nr, p in enumerate(paths) if p.name in tn][::-1]  # ids of tn paths, reversed
                if len(wtc.spheres):
                    n = max(1, optimal_threads.threads_count)
                    with clui.pbar(maxval=len(range(0, len(tnpaths), n)), mess="AutoBarber in action:") as pbar:
                        Reader.reset()
                        pool = Pool(processes=optimal_threads.threads_count)
                        bp = partial(barber_paths, spheres=wtc.spheres, only_for_names=tn)

                        new_paths.reinit(pbar)
                        nr = 0
                        while len(tnpaths):
                            pool.apply_async(bp, args=([paths.pop(nn) for nn in tnpaths[:n]],),
                                             callback=new_paths.callback_cric_next)
                            tnpaths = tnpaths[n:]
                            nr += 1
                            if nr % n == 0:
                                gc.collect()
                        pool.close()
                        pool.join()
                    gc.collect()
                else:
                    clui.message('AutoBarber procedure skip, no spheres detected.')
                    for nn in tnpaths:
                        new_paths.paths.append(paths.pop(nn))
                    tnpaths = []
        paths = new_paths.paths
        if Reader.edges:
            with clui.pbar(len(paths), 'Waterfall fall:') as pbar:
                waterfall_me(paths, pbar)

    ######################################################################
    # following procedures are run only if autobarber
    ######################################################################

    if options.auto_barber:

        clui.message("Recreate separate paths:")
        pbar = clui.pbar(len(paths))
        Reader.reset()
        pool = Pool(processes=optimal_threads.threads_count)
        # ysp = partial(yield_single_paths, progress=True, passing=options.allow_passing_paths)
        n = int(max(1, len(paths) / optimal_threads.threads_count / 3))
        # n = max(1, optimal_threads.threads_count)
        spaths = []
        nr_all = 0
        for sps_nrs in map(ysp, (paths[i:i + n] for i in range(0, len(paths), n))):
            # for sps_nrs in pool.imap_unordered(ysp, (paths[i:i + n] for i in xrange(0, len(paths), n))):
            for sp, nr in sps_nrs:
                spaths.append(sp)
                pbar.update(nr_all + nr + 1)
            nr_all += nr + 1
        pool.close()
        pool.join()
        gc.collect()
        pbar.finish()

        if Reader.edges:
            with clui.pbar(len(spaths), "Clean IDs:") as pbar:
                correct_spaths_ids(spaths, pbar)

        clui.message("(Re)Created %d separate paths out of %d raw paths" % (len(spaths), len(paths)))

    ######################################################################

    if options.auto_barber:

        if options.discard_short_paths or options.discard_short_object:
            if short_paths is not None or short_object is not None:
                with clui.pbar(len(spaths), discard_message) as pbar:
                    spaths_nr = len(spaths)
                    Reader.reset()
                    pool = Pool(processes=optimal_threads.threads_count)
                    # dse = partial(discard_short_etc, short_paths=short_paths, short_object=short_object,
                    #              short_logic=short_logic)
                    n = max(1, optimal_threads.threads_count)
                    new_spaths = NP(None)
                    new_spaths.reinit(pbar)
                    # CRIC AWARE MP!
                    if short_object is not None:
                        Reader.reset()
                        for pid in list(spaths_paths_ids(spaths)):
                            pool.apply_async(dse, args=(
                                [spaths.pop(nr) for nr in range(len(spaths) - 1, -1, -1) if spaths[nr].id.id == pid],),
                                             callback=new_spaths.callback_cric_next)
                    else:
                        for pid in list(spaths_paths_ids(spaths)):
                            pool.apply_async(dse, args=(
                                [spaths.pop(nr) for nr in range(len(spaths) - 1, -1, -1) if spaths[nr].id.id == pid],),
                                             callback=new_spaths.callback_next)
                    pool.close()
                    pool.join()
                    save_cric()
                    gc.collect()
                    spaths = new_spaths.paths
                    del new_spaths
                # del spaths
                spaths_nr_new = len(spaths)
                if spaths_nr == spaths_nr_new:
                    clui.message("No paths were discarded.")
                else:
                    clui.message("%d paths were discarded." % (spaths_nr - spaths_nr_new))
                    if options.remove_unused_parts:
                        with clui.pbar(len(spaths) + len(paths), "Removing (again) unused parts of paths:") as pbar:
                            paths = yield_generic_paths(spaths, progress=pbar)
                        if Reader.edges:
                            with clui.pbar(len(paths), 'Waterfall fall:') as pbar:
                                waterfall_me(paths, pbar)
            else:
                clui.message("No paths were discarded - no values were set.")

    ######################################################################

    if options.sort_by_id:
        with clui.fbm("Sort separate paths by resid"):
            spaths = sorted(spaths, key=lambda sp: (sp.id.id, sp.id.nr))

    ######################################################################

    # center of object
    coos = None
    if options.calculate_coo:
        with clui.pbar(len(spaths), "Center of object calculation") as pbar:
            spaths_nr = len(spaths)
            Reader.reset()
            pool = Pool(processes=optimal_threads.threads_count)

            n = max(1, optimal_threads.threads_count)

            # coos = pool.imap_unordered(center_of_object, (spaths[i:i + n] for i in xrange(0, len(spaths), n)))
            # coos = imap(center_of_object, (sp for sp in spaths if not isinstance(sp,PassingPath)))
            coos = pool.imap_unordered(center_of_object, (sp for sp in spaths if not isinstance(sp, PassingPath)))

            # CRIC AWARE MP!
            Reader.reset()
            coos = [coo for nr, coo, cric in coos if (pbar.next(step=nr) is None) and (CRIC.update_cric(cric) is None)]
            # coos = list(chain.from_iterable((coo for nr, coo, cric in coos if
            #                                 (pbar.next(step=nr) is None) and (
            #                                         CRIC.update_cric(cric) is None))))
            save_cric()
            pool.close()
            pool.join()
            gc.collect()

            coos = np.array(coos).mean(0)

            logger.info('Center of object is %0.2f, %0.2f, %0.2f' % tuple(coos))

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

    # return {'paths': paths, 'spaths': spaths, 'options': options._asdict(), 'soptions': soptions._asdict()}
    return {'paths': paths, 'spaths': spaths, 'options': options._asdict(), 'soptions': soptions._asdict(),
            'center_of_object': coos}


################################################################################


# inlets_clustering
def stage_IV_run(config, options,
                 spaths=None,
                 center_of_system=None,
                 center_of_object=None,
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
    inlets_center = center_of_system
    if options.inlets_center == 'coo' and center_of_object is not None:
        inlets_center = center_of_object
    elif options.inlets_center != 'cos':
        logger.warning("Unknown value of `inlets_center` option %r. Falling back to 'cos'." % options.inlets_center)
    if options.inlets_center == 'coo' and center_of_object is None:
        logger.warning(
            "Value of 'coo' not available. Use `calculate_coo` in stage III. Falling back to 'cos'." % inlets_center)
    inls = Inlets(spaths, center_of_system=inlets_center,
                  passing=not options.exclude_passing_in_clustering,
                  pbar=pbar)
    pbar.finish()
    clui.message("Number of inlets: %d" % inls.size)

    def noo():
        # returns number of outliers
        if 0 in inls.clusters_list:
            return inls.clusters.count(0)
        return 0

    def clustering_clustering():
        # ***** CLUSTERING *****
        with clui.fbm("Performing clustering", cont=False):
            potentially_recursive_clustering(config=config,
                                             clustering_name=config.cluster_name(),
                                             inlets_object=inls,
                                             spaths=spaths,
                                             message='clustering',
                                             max_level=max_level)
        clui.message('Number of outliers: %d' % noo())

    def outliers_detection():
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

    def reclustering():
        # ***** RECLUSTERING *****
        if options.recluster_outliers:
            with clui.fbm("Performing reclustering of outliers", cont=False):
                clustering_function = get_clustering_method(rcoptions, config)
                # perform reclustering
                # inls.recluster_outliers(clustering_function)
                inls.recluster_cluster(clustering_function, 0)
            clui.message('Number of clusters detected so far: %d' % len(inls.clusters_list))
            clui.message('Number of outliers: %d' % noo())

    def singletons_removal():
        # ***** SINGLETONS REMOVAL *****
        if options.singletons_outliers:
            if int(options.singletons_outliers):
                with clui.fbm("Removing clusters of size %d" % int(options.singletons_outliers)):
                    inls.small_clusters_to_outliers(int(options.singletons_outliers))
                clui.message('Number of clusters detected so far: %d' % len(inls.clusters_list))
                clui.message('Number of outliers: %d' % noo())

        # TODO: Move it after master paths!

    def add_passing_paths_to_clusters():
        # ***** ADD PASSING PATHS TO CLUSTERS *****
        if options.exclude_passing_in_clustering and options.add_passing_to_clusters:
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
                spaths_passing_ids = [nr for nr in range(len(spaths)) if spaths[nr].is_passing()]
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

    def join_clusters():
        # ***** JOIN CLUSTERS *****
        if options.join_clusters:
            with clui.fbm("Join clusters") as emess:
                for c2j in options.join_clusters.split():
                    emess('%s' % c2j)
                    c2j = list(map(int, c2j.split('+')))
                    inls.join_clusters(c2j)

    def renumber_clusters():
        # ***** RENUMBER CLUSTERS *****
        if options.renumber_clusters:
            with clui.fbm("Renumber clusters"):
                inls.renumber_clusters()

    def remove_inlets_in_specified_clusters():
        # ***** REMOVE INLETS IN SPECIFIED CLUSTERS *****
        if options.remove_inlets:
            with clui.fbm("Remove inlets") as emess:
                inlets_removed_reference = []
                inlets_removed_type = []
                for i2r in map(int, options.remove_inlets.split()):  # i2r is cluster id
                    emess('%d' % i2r)
                    inlets_removed_reference.extend(inls.lim2clusters(i2r).refs)
                    inlets_removed_type.extend(inls.lim2clusters(i2r).types)
                    inls.remove_inlets(i2r)
                # clusters removed, modify spaths
                for r, t in zip(inlets_removed_reference, inlets_removed_type):
                    # find path of r
                    p = next((sp for sp in spaths if sp.id == r))
                    p.remove_inlet(t)

    if inls.size > 0:
        clustering_order = [clustering_clustering, reclustering, singletons_removal, add_passing_paths_to_clusters,
                            join_clusters, renumber_clusters, remove_inlets_in_specified_clusters]
        if options.clustering_order == "aquarius":
            clui.message("Age of Aquarius.")
            clustering_order = [clustering_clustering, join_clusters, renumber_clusters, outliers_detection,
                                reclustering, singletons_removal, add_passing_paths_to_clusters,
                                remove_inlets_in_specified_clusters]
        elif options.clustering_order != "old-school":
            logger.waring("Unknown clustering order mode '%s', falling back to 'old-school'")

        for cluster_function in clustering_order:
            cluster_function()
            gc.collect()

        ###################
        # CLUSTERING DONE #
        ###################

        # ***** CLUSTERS' HISTORY *****
        clui.message('Clustering history:')
        clui.message(clui.print_simple_tree(inls.tree, prefix='').rstrip())

        traced_names = get_traced_names(spaths)

        def iter_over_tn():
            if options.separate_master:
                for tn in traced_names:
                    yield [tn]
                if options.separate_master_all:
                    yield traced_names
            else:
                yield traced_names

        # but only if user wants this
        master_paths = {}  # this and following dict hold master paths, keys here ar ctypes - one ctype one master path
        master_paths_smooth = {}
        if options.create_master_paths:
            fof = lambda sp: make_fractionof(sp, f=options.master_paths_amount)
            for tn in iter_over_tn():  # tn is a list of traced names
                # here we can select which paths are to be used
                mpsps = [sp for sp in spaths if
                         not isinstance(sp, PassingPath) and sp.id.name in tn]  # no PassingPaths, tn only
                # sort by size to ger stratification like effect
                mpsps = sorted(fof(mpsps), key=lambda sp: sp.size)
                ctypes = inls.spaths2ctypes(mpsps)  # temp ctypes
                with clui.fbm("Master paths calculations for %s" % (', '.join(tn)), cont=False):
                    smooth = get_smooth_method(soptions)  # this have to preceed GCS
                    if GCS.cachedir or GCS.cachemem:
                        pbar = clui.pbar(len(mpsps) * 2, mess='Building coords cache')  # +2 to be on the safe side
                        # TODO: do it in parallel
                        [sp.get_coords(smooth=None) for sp in mpsps if
                         pbar.next() is None and not isinstance(sp, PassingPath)]
                        [sp.get_coords(smooth=smooth) for sp in mpsps if
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

                        pbar = clui.pbar(len(mpsps) * 2)
                        for nr, ct in enumerate(ctypes_generic_list):
                            logger.debug('CType %s (%d)' % (str(ct), nr))
                            sps = lind(mpsps, what2what(ctypes_generic, [ct]))
                            if not len(sps):
                                logger.debug(
                                    'CType %s (%d), no single paths found, MasterPath calculation skipped.' % (
                                        str(ct), nr,))
                                continue
                            logger.debug('CType %s (%d), number of spaths %d' % (str(ct), nr, len(sps)))
                            ctspc = CTypeSpathsCollection(spaths=sps, ctype=ct, pbar=pbar,
                                                          threads=use_threads)
                            if tn == traced_names:  # either no separate or all option
                                master_paths.update({ct: ctspc.get_master_path(resid=(0, nr))})
                                master_paths_smooth.update({ct: ctspc.get_master_path(resid=(0, nr), smooth=smooth)})
                            elif len(tn) == 1:  # tn is one name
                                if tn[0] not in master_paths:
                                    master_paths.update({tn[0]: {}})
                                if tn[0] not in master_paths_smooth:
                                    master_paths_smooth.update({tn[0]: {}})
                                master_paths[tn[0]].update({ct: ctspc.get_master_path(resid=(0, nr))})
                                master_paths_smooth[tn[0]].update(
                                    {ct: ctspc.get_master_path(resid=(0, nr), smooth=smooth)})
                            else:
                                raise AssertionError(
                                    "Internal bug, not consistent iteration over traced names. Please send bug report to <%s>" % __mail__)
                            del ctspc
                        pbar.finish()
        gc.collect()
    else:
        clui.message("No inlets found. Clustering skipped.")
        # make empty results
        master_paths = {}
        master_paths_smooth = {}

    # ***** CLUSTERS' TYPES *****
    with clui.fbm("Calculating cluster types"):
        ctypes = inls.spaths2ctypes(spaths)

    ################################################################################

    return {'inls': inls,
            'ctypes': ctypes,
            'master_paths': master_paths,
            'master_paths_smooth': master_paths_smooth,
            'spaths': spaths}


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
    pa('Aqua-Duct: v%s' % aquaduct_version_nice())
    pa('NumPy: v%s' % np.__version__)
    pa('MDAnalysis: v%s' % mda_ver())

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
    traced_names = get_traced_names(spaths)

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
        yield traced_names, ''  # all types
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

        pbar = clui.pbar(maxval=nr_tnspt * 3 + nr_f * 2 * len(spaths) + len(spaths), mess="Calculating stats:")

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

        if Reader.sandwich_mode:
            for layer in range(Reader.number_of_layers()):
                pa("Clusters summary - inlets%s for layer %d" % (message, layer))
                pa.thead(header_line)
                for nr, cl in enumerate(inls.clusters_list):
                    inls_lim = inls.lim2spaths(
                        [sp for sp in spaths if isinstance(sp, sptype) and sp.id.id[0] == layer]). \
                        lim2rnames(tname). \
                        lim2clusters(cl)
                    pa(make_line(line_template, clusters_inlets(cl, inls_lim)), nr=nr)
                pa.tend(header_line)

        if pbar:
            pbar.next()

    ############
    if options.cluster_area:
        pa.sep()
        header_line, line_template = get_header_line_and_line_template(clusters_area_header(), head_nr=head_nr)
        for tname, sptype, message in iter_over_tnspt():
            pa("Clusters summary - areas%s" % message)
            pa.thead(header_line)
            for nr, cl in enumerate(inls.clusters_list):

                inls_lim = inls.lim2spaths([sp for sp in spaths if isinstance(sp, sptype)]).lim2rnames(
                    tname).lim2clusters(
                    cl)
                if inls_lim.size < 3: continue
                pa(make_line(line_template, clusters_area(cl, inls_lim, points=float(options.cluster_area_precision),
                                                          expand_by=float(options.cluster_area_expand))), nr=nr)
            pa.tend(header_line)

            if Reader.sandwich_mode:
                for layer in range(Reader.number_of_layers()):
                    pa("Clusters summary - areas%s for layer %d" % (message, layer))
                    pa.thead(header_line)
                    for nr, cl in enumerate(inls.clusters_list):

                        inls_lim = inls.lim2spaths(
                            [sp for sp in spaths if isinstance(sp, sptype) and sp.id.id[0] == layer]). \
                            lim2rnames(tname). \
                            lim2clusters(cl)
                        if inls_lim.size < 3: continue
                        pa(make_line(line_template,
                                     clusters_area(cl, inls_lim, points=float(options.cluster_area_precision),
                                                   expand_by=float(options.cluster_area_expand))), nr=nr)

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

        if Reader.sandwich_mode:
            for layer in range(Reader.number_of_layers()):
                header_line, line_template = get_header_line_and_line_template(clusters_stats_prob_header(),
                                                                               head_nr=head_nr)
                pa("Clusters statistics (of paths%s) probabilities of transfers for layer %d" % (message, layer))
                pa.thead(header_line)
                for nr, cl in enumerate(inls.clusters_list):
                    sp_ct_lim = ((sp, ct) for sp, ct in zip(spaths, ctypes) if
                                 cl in ct.clusters and isinstance(sp, sptype) and sp.id.name in tname)
                    pa(make_line(line_template, clusters_stats_prob(cl, sp_ct_lim)), nr=nr)
                    if pbar:
                        pbar.heartbeat()
                pa.tend(header_line)

                header_line, line_template = get_header_line_and_line_template(clusters_stats_len_header(),
                                                                               head_nr=head_nr)
                pa("Clusters statistics (of paths%s) mean lengths of transfers for layer %d" % (message, layer))
                pa.thead(header_line)
                for nr, cl in enumerate(inls.clusters_list):
                    sp_ct_lim = ((sp, ct) for sp, ct in zip(spaths, ctypes) if
                                 cl in ct.clusters and isinstance(sp, sptype) and sp.id.name in tname)
                    pa(make_line(line_template, clusters_stats_len(cl, sp_ct_lim)), nr=nr)
                    if pbar:
                        pbar.heartbeat()
                pa.tend(header_line)

                header_line, line_template = get_header_line_and_line_template(clusters_stats_steps_header(),
                                                                               head_nr=head_nr)
                pa("Clusters statistics (of paths%s) mean frames numbers of transfers for layer %d" % (message, layer))
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

        if Reader.sandwich_mode:
            for layer in range(Reader.number_of_layers()):
                pa("Separate paths clusters types summary - mean lengths of paths%s for layer %d" % (message, layer))
                pa.thead(header_line)

                total_size = len([sp for sp in spaths if sp.id.name in tname and isinstance(sp, sptype)])
                for nr, ct in enumerate(ctypes_generic_list):
                    sps = lind(spaths, what2what(ctypes_generic, [ct]))
                    sps = [sp for sp in sps if sp.id.name in tname and isinstance(sp, sptype)]
                    # ctypes_size.append(len(sps))
                    if len(sps) > 0:
                        pa(make_line(line_template, ctypes_spaths_info(ct, sps, add_size_p100=total_size, show="len")),
                           nr=nr)
                    if pbar:
                        pbar.next(len(sps))
                pa.tend(header_line)
                pa("Separate paths clusters types summary - mean number of frames of paths%s for layer %d" % (
                    message, layer))
                pa.thead(header_line)
                for nr, ct in enumerate(ctypes_generic_list):
                    sps = lind(spaths, what2what(ctypes_generic, [ct]))
                    sps = [sp for sp in sps if sp.id.name in tname and isinstance(sp, sptype)]
                    # ctypes_size.append(len(sps))
                    if len(sps) > 0:
                        pa(make_line(line_template,
                                     ctypes_spaths_info(ct, sps, add_size_p100=total_size, show="frames")),
                           nr=nr)
                    if pbar:
                        pbar.next(len(sps))

        pa.tend(header_line)

    ############
    pa.sep()
    pa("List of separate paths and properties")
    header_line, line_template = get_header_line_and_line_template(spath_full_info_header(total=True), head_nr=head_nr)
    pa.thead(header_line)
    for nr, (sp, ctype) in enumerate(zip_longest(spaths, ctypes, fillvalue=None)):
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
        # this fuction is redefined here because of changed messages
        yield traced_names, 'amol'  # all types
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
                header += ['%s_%d' % (s, number) for s in ['scope_area', 'scope_volume', 'object_area', 'object_volume']]
                fmt += ['%0.3f', '%0.2f'] * 2
            elif not len(scope_size):
                header += ['%s' % s for s in ['scope_area', 'scope_volume', 'object_area', 'object_volume']]
                fmt += ['%0.3f', '%0.2f'] * 2
            if Reader.sandwich_mode or not len(scope_size):
                scope_size.append([])
                object_size.append([])
            for frame in traj_reader.iterate_over_frames():
                scope = traj_reader.parse_selection(options.scope_chull)
                ch = scope.chull(inflate=options.scope_chull_inflate)
                if ch is not None:
                    scope_size[-1].append((ch.area, ch.volume))
                else:
                    logger.warning("Cannot get Convex Hull for Scope in %d:%d" % (number, frame))
                    scope_size[-1].append((float('nan'), float('nan')))
                res = traj_reader.parse_selection(options.object_chull)
                ch = res.chull()
                if ch is not None:
                    object_size[-1].append((ch.area, ch.volume))
                else:
                    logger.warning("Cannot get Convex Hull for Object in %d:%d" % (number, frame))
                    object_size[-1].append((float('nan'), float('nan')))
                pbar.next()
        for s_s, o_s in zip(scope_size, object_size):
            h = np.hstack((h, s_s, o_s))
        pbar.finish()
    # add frame column?
    frame_col = np.array([list(range(max_frame))]).T
    h = np.hstack((frame_col, h))
    header = ['frame'] + header
    fmt = ['%u'] + fmt
    # save???
    if options.save:
        h_fname = options.save + '.csv'
    else:
        import io as StringIO
        h_fname = StringIO.StringIO()
    np.savetxt(h_fname, h,
               fmt=fmt,
               delimiter=',',
               header=','.join(header))
    if not options.save:
        h_fname.getvalue()
        # print h_fname.getvalue()

    return {'hist': h, 'header': header}


################################################################################

# visualize


def stage_VI_run(config, options,
                 spaths=None,
                 inls=None,
                 ctypes=None,
                 master_paths=None,
                 master_paths_smooth=None,
                 center_of_system=None,
                 center_of_object=None,
                 **kwargs):
    # enable real cache of ort

    Reader.reset()

    from aquaduct.visual.pymol_connector import ConnectToPymol, SinglePathPlotter
    # from aquaduct.visual.pymol_connector import cmd as pymol_cmd
    from aquaduct.visual.helpers import ColorMapDistMap

    soptions = config.get_smooth_options()
    smooth = get_smooth_method(soptions)

    traced_names = get_traced_names(spaths)

    def iter_over_tn():
        if not options.split_by_type:
            yield None, ""
        else:
            if options.retain_all_types:
                yield None, ""
            for tn in traced_names:
                yield tn, "_%s" % str(tn)

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
                    if chull is not None:
                        spp.convexhull(chull, name='scope_shape%d' % nr, state=frame + 1)
                    else:
                        logger.debug('Scope convex hull calculations failed for frame %d.' % frame)

    if options.show_object_chull:
        with clui.fbm("Object shape"):
            for nr, traj_reader in enumerate(Reader.iterate(threads=False)):
                traj_reader = traj_reader.open()
                frames_to_show = range2int(options.show_object_chull_frames)
                for frame in frames_to_show:
                    traj_reader.set_frame(frame)
                    object_shape = traj_reader.parse_selection(options.show_object_chull)
                    chull = object_shape.chull()
                    if chull is not None:
                        spp.convexhull(chull, name='object_shape%d' % nr, color=np.array([255, 153, 0]) / 255.,
                                       state=frame + 1)  # orange
                    else:
                        logger.debug('Object convex hull calculations failed for frame %d.' % frame)

    fof = lambda sp: np.array(list(make_fractionof(sp, f=options.inlets_clusters_amount)))
    if options.inlets_clusters:
        tn_lim = lambda tn: inls if tn is None else inls.lim2rnames(tn)
        with clui.fbm("Clusters"):
            for tn, tn_name in iter_over_tn():
                # TODO: require stage V for that?
                # no_of_clusters = len(inls.clusters_list)  # total, including outliers
                cmap = ColorMapDistMap()
                for c in tn_lim(tn).clusters_list:
                    if Reader.sandwich_mode:
                        for layer in range(Reader.number_of_layers()):
                            # coords for current cluster and layer
                            ics = tn_lim(tn).lim_to([spath.id.id[0] for spath in spaths], layer).lim2clusters(c).coords
                            # Skip empty clusters
                            if not len(ics):
                                continue
                            if c == 0:
                                c_name = 'out'
                            else:
                                c_name = str(int(c))
                            spp.scatter(fof(ics), color=cmap(c), name="cluster_%s%s_L%s" % (c_name, tn_name, layer))
                            if False:  # TODO: This does not work any more in that way. Rewrite it or remove it
                                radii = tn_lim(tn).lim2clusters(c).radii
                                if len(radii) > 0:
                                    spp.scatter(ics, color=cmap(c), radius=radii,
                                                name="cluster_radii_%s%s" % (c_name, tn_name))
                    else:
                        # coords for current cluster
                        ics = tn_lim(tn).lim2clusters(c).coords
                        if c == 0:
                            c_name = 'out'
                        else:
                            c_name = str(int(c))
                        spp.scatter(fof(ics), color=cmap(c), name="cluster_%s%s" % (c_name, tn_name))
                        if False:  # TODO: This does not work any more in that way. Rewrite it or remove it
                            radii = tn_lim(tn).lim2clusters(c).radii
                            if len(radii) > 0:
                                spp.scatter(ics, color=cmap(c), radius=radii,
                                            name="cluster_radii_%s%s" % (c_name, tn_name))

    if options.cluster_area:
        from aquaduct.geom import hdr
        from aquaduct.geom.hdr_contour import hdr2contour, iscontour
        if iscontour:
            cmap = cmaps._cmap_jet_256
            with clui.fbm("Clusters contours"):
                for c in inls.clusters_list:
                    if Reader.sandwich_mode:
                        for layer in range(Reader.number_of_layers()):
                            # coords for current cluster
                            ics = inls.lim_to([spath.id.id[0] for spath in spaths], layer).lim2clusters(c).coords
                            if c == 0:
                                c_name = 'out'
                            else:
                                c_name = str(int(c))
                            # calcualte hdr
                            # print inls.center_of_system, c_name, len(ics), alt_center_of_system
                            if len(ics) < 3: continue
                            h = hdr.HDR(np.array(ics), points=float(options.cluster_area_precision),
                                        expand_by=float(options.cluster_area_expand),
                                        center_of_system=inls.center_of_system)
                            spp.multiline_begin()
                            for fraction in range(100, 0, -5):  # range(100, 85, -5) + range(80, 40, -10):
                                # print c_name + '_D%d' % fraction
                                coords = hdr2contour(h, fraction=fraction / 100.)
                                if coords is not None:
                                    color = cmap[int(255 * (1 - fraction / 100.))]
                                    spp.multiline_add(coords, color=color)
                            spp.multiline_end(name=c_name + '_L{}_DC'.format(layer))
                    else:
                        # coords for current cluster
                        ics = inls.lim2clusters(c).coords
                        if c == 0:
                            c_name = 'out'
                        else:
                            c_name = str(int(c))
                        # calcualte hdr
                        # print inls.center_of_system, c_name, len(ics), alt_center_of_system
                        if len(ics) < 3: continue
                        h = hdr.HDR(np.array(ics), points=float(options.cluster_area_precision),
                                    expand_by=float(options.cluster_area_expand),
                                    center_of_system=inls.center_of_system)
                        spp.multiline_begin()
                        for fraction in range(100, 0, -5):  # range(100, 85, -5) + range(80, 40, -10):
                            # print c_name + '_D%d' % fraction
                            coords = hdr2contour(h, fraction=fraction / 100.)
                            if coords is not None:
                                color = cmap[int(255 * (1 - fraction / 100.))]
                                spp.multiline_add(coords, color=color)
                        spp.multiline_end(name=c_name + '_DC')
            spp.scatter(np.array([center_of_system]), color=cmap[10], name="CoS")
            if center_of_object is not None:
                spp.scatter(np.array([center_of_object]), color=cmap[10], name="CoO")

    fof = lambda sp: list(make_fractionof(sp, f=options.ctypes_amount))

    # master paths can have some keys which are not of ct type - names of molecules
    # print master_paths.keys()
    # print master_paths_smooth.keys()
    master_paths_separate = [k for k in master_paths.keys() if isinstance(k, str)]

    # TODO: is isinstance good in this instance?
    if options.ctypes_raw:
        with clui.fbm("CTypes raw"):
            for nr, ct in enumerate(ctypes_generic_list):
                clui.message(str(ct), cont=True)
                for tn, tn_name in iter_over_tn():
                    if Reader.sandwich_mode:
                        for layer in range(Reader.number_of_layers()):
                            clui.message("Paths in layer {}/{}:".format(layer, Reader.number_of_layers() - 1))
                            sps = lind(spaths, what2what(ctypes_generic, [ct]))
                            tn_lim = lambda tn: sps if tn is None else [sp for sp in sps if tn == sp.id.name]
                            plot_spaths_traces(fof([spath for spath in tn_lim(tn) if spath.id.id[0] == layer]),
                                               name=str(ct) + '_raw' + tn_name + '_L' + str(layer),
                                               split=False,
                                               spp=spp)
                    else:
                        sps = lind(spaths, what2what(ctypes_generic, [ct]))
                        tn_lim = lambda tn: sps if tn is None else [sp for sp in sps if tn == sp.id.name]
                        plot_spaths_traces(fof(tn_lim(tn)), name=str(ct) + '_raw' + tn_name, split=False, spp=spp)
                for mp_nr in range(len(master_paths_separate) + 1):
                    mp_name = ""
                    mp = None
                    if mp_nr:
                        mp_name = master_paths_separate[mp_nr - 1]
                        if ct in master_paths[mp_name]:
                            mp = master_paths[mp_name][ct]
                            mp_name = "_" + mp_name
                    else:
                        if ct in master_paths:
                            mp = master_paths[ct]
                        else:
                            mp = None
                    if mp is not None:
                        plot_spaths_traces([mp],
                                           name=str(ct) + '_raw_master' + mp_name,
                                           split=False,
                                           spp=spp,
                                           smooth=lambda anything: anything)

    master_paths_separate = [k for k in master_paths_smooth.keys() if isinstance(k, str)]
    if options.ctypes_smooth:
        with clui.fbm("CTypes smooth"):
            for nr, ct in enumerate(ctypes_generic_list):
                clui.message(str(ct), cont=True)
                for tn, tn_name in iter_over_tn():
                    if Reader.sandwich_mode:
                        for layer in range(Reader.number_of_layers()):
                            clui.message("Paths in layer {}/{}:".format(layer, Reader.number_of_layers() - 1))

                            sps = lind(spaths, what2what(ctypes_generic, [ct]))
                            tn_lim = lambda tn: sps if tn is None else [sp for sp in sps if tn == sp.id.name]

                            plot_spaths_traces(fof([spath for spath in tn_lim(tn) if spath.id.id[0] == layer]),
                                               name=str(ct) + '_smooth' + tn_name + '_L' + str(layer),
                                               split=False,
                                               spp=spp,
                                               smooth=smooth)
                    else:
                        sps = lind(spaths, what2what(ctypes_generic, [ct]))
                        tn_lim = lambda tn: sps if tn is None else [sp for sp in sps if tn == sp.id.name]
                        plot_spaths_traces(fof(tn_lim(tn)), name=str(ct) + '_smooth' + tn_name, split=False, spp=spp,
                                           smooth=smooth)
                for mp_nr in range(len(master_paths_separate) + 1):
                    mp_name = ""
                    mp = None
                    if mp_nr:
                        mp_name = master_paths_separate[mp_nr - 1]
                        if ct in master_paths_smooth[mp_name]:
                            mp = master_paths_smooth[mp_name][ct]
                            mp_name = "_" + mp_name
                    else:
                        if ct in master_paths_smooth:
                            mp = master_paths_smooth[ct]
                        else:
                            mp = None
                    if mp is not None:
                        plot_spaths_traces([mp],
                                           name=str(ct) + '_smooth_master' + mp_name,
                                           split=False,
                                           spp=spp,
                                           smooth=lambda anything: anything)
                    mp = None
                    if mp_nr:
                        mp_name = master_paths_separate[mp_nr - 1]
                        if ct in master_paths[mp_name]:
                            mp = master_paths[mp_name][ct]
                            mp_name = "_" + mp_name
                    else:
                        if ct in master_paths_smooth:
                            mp = master_paths_smooth[ct]
                        else:
                            mp = None
                    if mp is not None:
                        plot_spaths_traces([mp],
                                           name=str(ct) + '_raw_master_smooth' + mp_name,
                                           split=False,
                                           spp=spp,
                                           smooth=smooth)

    def plot_paths(paths, tn_name, layer_name=""):
        if options.all_paths_raw:
            with clui.fbm("All raw paths" + tn_name.replace('_', ' ')):
                plot_spaths_traces(fof(paths), name='all_raw' + tn_name + layer_name,
                                   split=options.all_paths_split, spp=spp)
        if options.all_paths_raw_io:
            with clui.fbm("All raw paths io" + tn_name.replace('_', ' ')):
                plot_spaths_inlets(fof(paths), name='all_raw_paths_io' + tn_name + layer_name, spp=spp)

        if options.all_paths_smooth:
            with clui.fbm("All smooth paths" + tn_name.replace('_', ' ')):
                plot_spaths_traces(fof(paths), name='all_smooth' + tn_name + layer_name,
                                   split=options.all_paths_split, spp=spp,
                                   smooth=smooth)
        if options.all_paths_smooth_io:
            with clui.fbm("All smooth paths io" + tn_name.replace('_', ' ')):
                plot_spaths_inlets(fof(paths), name='all_smooth_paths_io' + tn_name + layer_name, spp=spp)

        with clui.fbm("Paths as states" + tn_name.replace('_', ' ')):
            if options.paths_raw:
                clui.message("raw", cont=True)
                plot_spaths_traces(paths, name='raw_paths' + tn_name + layer_name,
                                   states=options.paths_states,
                                   separate=not options.paths_states,
                                   spp=spp)
            if options.paths_smooth:
                clui.message("smooth", cont=True)
                plot_spaths_traces(paths, name='smooth_paths' + tn_name + layer_name,
                                   states=options.paths_states,
                                   separate=not options.paths_states, smooth=smooth, spp=spp)
            if options.paths_raw_io:
                clui.message("raw_io", cont=True)
                plot_spaths_inlets(paths, name='raw_paths_io' + tn_name + layer_name,
                                   states=options.paths_states,
                                   separate=not options.paths_states, spp=spp)
            if options.paths_smooth_io:
                clui.message("smooth_io", cont=True)
                plot_spaths_inlets(paths, name='smooth_paths_io' + tn_name + layer_name,
                                   states=options.paths_states,
                                   separate=not options.paths_states, smooth=smooth, spp=spp)

    fof = lambda sp: list(make_fractionof(sp, f=options.all_paths_amount))
    tn_lim = lambda tn: spaths if tn is None else [sp for sp in spaths if tn == sp.id.name]

    for tn, tn_name in iter_over_tn():
        if Reader.sandwich_mode:
            for layer in range(Reader.number_of_layers()):
                clui.message("Paths in layer {}/{}:".format(layer, Reader.number_of_layers() - 1))
                plot_paths([spath for spath in tn_lim(tn) if spath.id.id[0] == layer], tn_name, "_L{}".format(layer))
        else:
            plot_paths(tn_lim(tn), tn_name)

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


__all__ = '''valve_exec_stage
stage_I_run
stage_II_run
stage_III_run
stage_IV_run
stage_V_run
stage_VI_run
'''.split()
