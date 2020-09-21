#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2018-2019  Tomasz Magdziarz
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

"""
What have I got in my pocket?
"""



################################################################################
# reuse AQ logger

from __future__ import absolute_import
from __future__ import print_function
import logging
from aquaduct import logger, logger_name

formatter_string = '%(name)s:%(levelname)s:[%(module)s|%(funcName)s@%(lineno)d]: %(message)s'
# create and add console handler with WARNING level to the AQ logger
formatter = logging.Formatter(formatter_string)
ch = logging.StreamHandler()
ch.setFormatter(formatter)
ch.setLevel(logging.WARNING)  # default level is WARNING
logger.addHandler(ch)

################################################################################

import json
import gzip
import csv
import sys
from collections import namedtuple

################################################################################

# Boltzmann k constant
# k = 0.0019872041 # kcal/mol/K
k = 0.0083144621  # kJ/mol/K
k_unit = 'kJ/mol/K'

################################################################################

class count_ref(object):

    def __init__(self, pbar, ref_sel):
        self.pbar = pbar  # queue
        self.ref_sel = ref_sel

    def __call__(self, traj_reader):
        traj_reader = traj_reader.open()
        ref = 0.
        for frame in traj_reader.iterate():
            ref += len(list(traj_reader.parse_selection(self.ref_sel).residues().ids()))
            self.pbar.put(1)  # report one frame a time
        return ref


################################################################################


if __name__ == "__main__":

    from sys import exc_info

    try:
        from aquaduct.utils import clui

        with clui.tictoc('What have I got in my pocket?'):

            # ----------------------------------------------------------------------#
            # argument parsing

            import argparse
            from aquaduct import version_nice as aquaduct_version_nice

            description_version = '''Aquaduct library version %s''' % (aquaduct_version_nice(),)
            description = '''What have I got in my pocket?'''

            parser = argparse.ArgumentParser(description=description,
                                             formatter_class=argparse.ArgumentDefaultsHelpFormatter)

            parser.add_argument("-c", action="store", dest="config_file", required=False,
                                help="Config file filename.")
            parser.add_argument("-t", action="store", dest="threads", required=False, default=None,
                                help="Limit Aqua-Duct calculations to given number of threads.")
            parser.add_argument("-r", action="store", dest="results_dir", required=False,
                                help="Path to results directory", default="", type=str)
            parser.add_argument("--debug", action="store_true", dest="debug", required=False,
                                help="Prints debug info.")
            parser.add_argument("--debug-file", action="store", dest="debug_file", required=False,
                                help="Debug log file.")
            parser.add_argument("--paths-types", action="store", dest="paths_types", type=str, required=False,
                                default="",
                                help="Limit calculations to given paths types, i.e. given molecules.")
            parser.add_argument("--raw", action="store_true", dest="raw", required=False,
                                help="Use raw data from paths instead of single paths.")
            parser.add_argument("--raw-master", action="store_true", dest="raw_master", required=False,
                                help="Use raw data from paths instead of single paths, only in master paths calculations.")
            parser.add_argument("--raw-discard-singletons", action="store", dest="raw_singl", required=False, type=int,
                                default=1,
                                help="Discard short scope only segments from raw data.")
            parser.add_argument("--window-full", action="store_true", dest="wfull", required=False,
                                help="Return full window if windows is used.")
            parser.add_argument("--windows", action="store", dest="windows", type=int, required=False, default=1,
                                help="Number of windows to calculate.")
            parser.add_argument("--wsize", action="store", dest="wsize", type=int, required=False, default=None,
                                help="Size of window in frames.")
            parser.add_argument("--reference-value", action="store", dest="ref", type=float, required=False, default=None,
                                help="Reference value in [%s]." % k_unit)
            parser.add_argument("--reference-calc", action="store_true", dest="ref_calc", required=False,
                                help="Calculate reference value with scope and reference molecules.")
            parser.add_argument("--reference-radius", action="store", dest="ref_radius", type=float, required=False,
                                default=2.,
                                help="Radius of reference in [Å].")
            parser.add_argument("--reference-mol", action="store", dest="ref_mol", type=str, required=False,
                                default='resname WAT',
                                help="Selection of reference molecules.")
            parser.add_argument("--temperature", action="store", dest="temp", type=float, required=False, default=300.,
                                help="Simulation temperature in [K].")
            parser.add_argument("--gsize", action="store", dest="grid_size", type=float, required=False, default=1.,
                                help="Size of grid's cells in in [Å].")
            parser.add_argument("--pockets", action="store_true", dest="pockets", required=False,
                                help="Calculate pockets.")
            parser.add_argument("--hotspots", action="store_true", dest="hotspots", required=False,
                                help="Calculates hotspots if pockets are calculated.")
            parser.add_argument("--energy-profile", action="store_true", dest="eng_pro", required=False,
                                help="Calculates energy profiles for master paths.")
            # FIXME: Temporary option for loading master paths when they are not needed due to master_radius default value.
            parser.add_argument("--master", action="store_true", dest="master", required=False,
                                help="Enables master paths calculation.")
            parser.add_argument("--master-radius", action="store", dest="master_radius", type=float, required=False, default=2.,
                                help="Calculate profiles for master paths with given radius in [Å].")
            parser.add_argument("--master-ctypes", action="store", dest="master_ctypes", type=str, required=False,
                                default="",
                                help="Limit calculations to given ctypes.")
            # parser.add_argument("--master-radii", action="store_true", dest="master_radii", required=False,
            #                    help="Calculate profiles for master paths using width as radii.")
            parser.add_argument("--io-threshold", action="store", dest="io_threshold", type=float, required=False,
                                default=None,
                                help="Percent value of maximal density which will be used to partition pocket into inner and outer instead of mean value.")
            parser.add_argument("--path-id", action="store", dest="path_id", type=str, required=False,
                                default=None, help="Calculate profiles for specified path ID.")
            parser.add_argument("--path-file", action="store", dest="path_file", type=str, required=False,
                                help="Use coordinates from specified CSV file.")
            parser.add_argument("--path-radius", action="store", dest="path_radius", type=float, required=False,
                                default=2., help="Calculate profiles for path with given radius in [Å].")
            parser.add_argument("--path-smooth", action="store_true", dest="path_smooth", required=False,
                                help="If used path coordinates will be smoothed.")
            parser.add_argument("--raw-path", action="store_true", dest="raw_path", required=False,
                                help="Use raw data from paths instead of single paths. "
                                     "Used for path energy profiles calculations and for extracting raw path.")
            parser.add_argument("--extract-path", action="store", dest="extract_path", type=str, required=False,
                                default=None, help="Extract path coordinates with specified ID.")
            parser.add_argument("--output-file", action="store", dest="output_file", type=str, required=False,
                                default="path_coords.csv", help="Output CSV filename for extracted coordinates.")
            parser.add_argument("--output-suffix", action="store", dest="output_suffix", type=str, required=False, default="", help="Arbitrary suffix appended to output filenames.")


            args = parser.parse_args()

            # ----------------------------------------------------------------------#

            # debug
            # at this stage logger is the AQ root logger
            if args.debug:
                logger.removeHandler(ch)  # remove old ch handlers
                CH = logging.StreamHandler()
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
            logger = logging.getLogger(logger_name + '.pond')
            logger.info('Initialization of Pond logging done.')

            # ----------------------------------------------------------------------#
            # config

            from aquaduct.apps.valve import ValveConfig, valve_load_config

            config = ValveConfig()  # config template
            valve_load_config(args.config_file, config)
            # get global options
            goptions = config.get_global_options()

            # ----------------------------------------------------------------------#
            # init and options

            from aquaduct.utils.multip import optimal_threads
            from aquaduct.traj.sandwich import Reader, Window
            from aquaduct.apps.valve.data import get_vda_reader
            from aquaduct.apps.data import save_cric
            from os import pathsep
            import numpy as np
            from multiprocessing import Pool, Manager
            from aquaduct.geom import pocket
            from aquaduct.traj.dumps import WriteMOL2
            from aquaduct.apps.valve.helpers import get_linearize_method
            from aquaduct.apps.valve.helpers import get_smooth_method
            from aquaduct.geom import traces
            
            import os
            from aquaduct.geom.smooth import SavgolSmooth

            if args.threads is None:
                optimal_threads.threads_count = optimal_threads.cpu_count + 1
            else:
                optimal_threads.threads_count = int(args.threads)
            clui.message("Number of threads Pond is allowed to use: %d" % optimal_threads.threads_count)
            if (1 < optimal_threads.threads_count < 3) or (
                    optimal_threads.threads_count - 1 > optimal_threads.cpu_count):
                clui.message(
                    "Number of threads is not optimal; CPU count reported by system: %d" % optimal_threads.cpu_count)
            # because it is used by mp.Pool it should be -1???
            if optimal_threads.threads_count > 1:
                optimal_threads.threads_count -= 1
                clui.message("Main process would use 1 thread.")
                clui.message("Concurent calculations would use %d threads." % optimal_threads.threads_count)

            # Maximal frame checks
            frames_window = Window(goptions.min_frame, goptions.max_frame, goptions.step_frame)

            # Open trajectory reader
            # trajectories are split with os.pathsep

            Reader(goptions.top, goptions.trj,
                   window=frames_window,
                   sandwich=goptions.sandwich,
                   threads=optimal_threads.threads_count)  # trajectory reader

            # ----------------------------------------------------------------------#
            # results dir
            rdir = ""
            if len(args.results_dir):
                if not os.path.exists(args.results_dir):
                    os.makedirs(args.results_dir)
                rdir = args.results_dir
            if len(rdir) == 0:
                rdir += '.'
            rdir += os.path.sep

            results_meta = {'options': vars(args)}
            rmu = lambda k, v: results_meta.update({k: v})
            args.output_suffix = "_%s" % args.output_suffix

            # ----------------------------------------------------------------------#
            # load paths
            paths_types = [pt.strip() for pt in args.paths_types.split(' ') if len(pt.strip())]
            ptn = ''  # additiona paths types name suffix
            if paths_types:
                rmu('paths_types', paths_types)
                clui.message('Limiting calculations to paths of %s.' % ', '.join(paths_types))
                ptn = '_' + '_'.join(paths_types)

            if (args.pockets or args.master_radius) and not (args.master_radius and not args.raw and args.raw_master):
                with clui.tictoc('Loading paths'):
                    if not args.raw:
                        # get stage III options
                        options3 = config.get_stage_options(2)
                        with clui.fbm('Loading data dump from %s file' % options3.dump):
                            vda = get_vda_reader(options3.dump, mode='r')
                            result3 = vda.load()
                        paths = result3.pop('spaths')
                        if paths_types:
                            paths = [p for p in paths if p.id.name in paths_types]
                        rmu('paths', 'spaths')
                    else:
                        # get stage II options
                        options2 = config.get_stage_options(1)
                        with clui.fbm('Loading data dump from %s file' % options2.dump):
                            vda = get_vda_reader(options2.dump, mode='r')
                            result2 = vda.load()
                        with clui.fbm('Preparing raw data'):
                            paths = result2.pop('paths')
                            if args.raw_singl:
                                for p in paths:
                                    p.discard_singletons(singl=args.raw_singl)
                            paths = [p for p in paths if len(p.frames)]
                            if paths_types:
                                paths = [p for p in paths if p.name in paths_types]
                        rmu('paths', 'paths')

            # ----------------------------------------------------------------------#
            # reference value

            ref = args.ref
            # calculate n_z for reference
            if args.ref_calc:
                with clui.tictoc('Reference calculation'):
                    from aquaduct.geom import traces
                    from scipy.spatial.distance import cdist
                    # 1. Get scope atoms.
                    # 2. Calculate convexhull.
                    # 3. Get first vertex and in the direction in which it is oriented look for farthest reference molecule.
                    # 4. Point for reference calculation get as midpoint between the used vertex and the farthest molecule.
                    # 5. Proceed normally.
                    options2 = config.get_stage_options(1)

                    for traj_reader in Reader.iterate():
                        traj_reader = traj_reader.open()
                        for frame in traj_reader.iterate():
                            with clui.pbar(mess="Automatic reference point calculation",maxval=7) as pbar:
                                # 1.
                                scope = traj_reader.parse_selection(options2.scope)
                                pbar.next()
                                # 2.
                                scope_ch = scope.chull(inflate=options2.scope_convexhull_inflate)
                                pbar.next()
                                # 3.
                                v = scope_ch.vertices_points[0] # vertex
                                c = np.mean(scope_ch.vertices_points) # center of scope
                                #ref_mol = traj_reader.parse_selection(args.ref_mol).residues().coords()
                                ref_mol = traj_reader.select_all().residues().coords()
                                ref_mol = np.array(list(ref_mol)) # reference molecules coordinates
                                pbar.next()
                                # find all ref_mol that are above the vertex
                                ref_mol = ref_mol[[traces.is_p_above_vp0_plane(rm,v-c,v) > 0 for rm in ref_mol]]
                                pbar.next()
                                # find distances of remaining ref_mol to vc line
                                rmd = [traces.distance_p_to_ab(rm,v,c) for rm in ref_mol]
                                pbar.next()
                                # find all ref_mol that are within reference radius from vc line
                                ref_mol = ref_mol[[d < args.ref_radius for d in rmd]]
                                pbar.next()
                                # find the most distant ref_mol
                                com = cdist(ref_mol,[v])
                                com = ref_mol[np.argmax(com)]
                                pbar.next()
                            # 4.
                            clui.message("Using reference point as middle between (%0.4f,%0.4f,%0.4f) and (%0.4f,%0.4f,%0.4f)." % (tuple(map(float,v))+tuple(map(float,com))))
                            com = (com+v)/2
                            clui.message("Using reference point (%0.4f,%0.4f,%0.4f)." % tuple(map(float,com)))
                            break
                        break
                    Reader.reset()

                    goal = Reader.number_of_frames(onelayer=False)
                    with clui.pbar(goal, mess="Calculating density in the reference area:") as pbar:
                        ref = []
                        ref_sel = '(%s) and (point %f %f %f %f)' % (
                                (args.ref_mol,) + tuple(map(float, com)) + (args.ref_radius,))
                        # ref_sel = '(%s) and (sphzone %f %s)' % (args.ref_mol, args.ref_radius, args.ref)
                        manager = Manager()
                        pbar_queue = manager.Queue()
                        pool = Pool(processes=optimal_threads.threads_count)

                        r = pool.map_async(
                            count_ref(pbar_queue, ref_sel),
                            Reader.iterate(),
                            callback=ref.extend)

                        progress = 0
                        for p in iter(pbar_queue.get, None):
                            progress += p
                            pbar.next(step=p)
                            if progress == goal:
                                break
                        r.wait()
                        pool.close()
                        pool.join()
                        del pbar_queue, manager

                    ref = float(sum(ref))

                    clui.message('Reference number of molecules: %d [molecules].' % int(ref))

                    ref /= float(Reader.number_of_frames(onelayer=False))  # if sandwich mean value will be correct

                    ref /= 4. / 3. * np.pi * (float(args.ref_radius) ** 3)
                    clui.message('Reference density: %0.4f [molecules/A^3].' % ref)

                    # Boltzmann k constant
                    # k = 0.0019872041 # kcal/mol/K
                    # k = 0.0083144621  # kJ/mol/K
                    rmu('energy_unit', k_unit)

                    ref = -k * args.temp * np.log(ref)
                    rmu('reference_correction', float(ref))
                    clui.message('Reference value: %0.4f [kJ/mol/K].' % ref)

            if ref:
                rmu('reference_density_correction', np.exp(ref/(-k * args.temp)))

            # ----------------------------------------------------------------------#
            # windows

            W = args.windows
            WS = args.wsize

            many_windows = W > 1 or (W == 1 and (WS is not None))
            many_windows = many_windows and WS < Reader.number_of_frames(onelayer=True)

            windows = []
            clui.message("Calculation will be done for following windows:")
            for wnr, window in enumerate(pocket.windows(Reader.number_of_frames(onelayer=True), windows=W, size=WS)):
                middle = sum(window) / 2
                if wnr and many_windows:
                    print(("W%d %d:%d middle: %d" % (wnr, window[0], window[1], middle)))
                    windows.append([wnr, window[0], window[1], middle])
                elif (wnr == 0) and (args.wfull or (not many_windows)):
                    print(("full %d:%d" % (window[0], window[1])))
                    rmu('full_window',[window[0], window[1]])
                rmu('windows', windows)

            # ----------------------------------------------------------------------#
            # calculate pockets

            if args.pockets:
                clui.message('What have I got in my pocket?')
                with clui.tictoc('Pockets & hotspots calculation'):

                    W = args.windows
                    WS = args.wsize
                    grid_size = args.grid_size
                    grid_area = grid_size ** 3
                    with clui.pbar(len(paths) * (1 + W + int(args.wfull)), mess='Calculating pockets:') as pbar:
                        pockets_volume = open(rdir + ('volumes%s%s.dat' % (ptn, args.output_suffix), 'w')
                        pockets_volume.write(('\t'.join('W_start W_end Outer Inner'.split())) + os.linesep)
                        pool = Pool(processes=optimal_threads.threads_count)
                        edges = pocket.find_edges(paths, grid_size=grid_size, pbar=pbar, map_fun=pool.imap_unordered)
                        number_of_frames = Reader.number_of_frames(onelayer=True)
                        if WS is None:
                            WSf = float(number_of_frames / float(W))
                        else:
                            WSf = float(WS)
                        if Reader.sandwich_mode:
                            WSf *= Reader.number_of_layers()

                        wmol2 = None
                        hsmol2 = None
                        for wnr, window in enumerate(
                                pocket.windows(Reader.number_of_frames(onelayer=True), windows=W, size=WS)):
                            number_of_frames = (window[-1] - window[0])
                            if Reader.sandwich_mode:
                                number_of_frames *= Reader.number_of_layers()

                            if wnr and many_windows:
                                if wmol2 is None:
                                    wmol2 = [WriteMOL2(rdir + 'outer%s%s.mol2' % (ptn, args.output_suffix)),
                                             WriteMOL2(rdir + 'inner%s%s.mol2' % (ptn, args.output_suffix))]
                                if hsmol2 is None and args.hotspots:
                                    hsmol2 = WriteMOL2(rdir + 'hotspots%s%s.mol2' % (ptn, args.output_suffix))
                                D = pocket.distribution(paths, grid_size=grid_size, edges=edges, window=window,
                                                        pbar=pbar, map_fun=pool.imap_unordered)
                                H = (D[-1] / WSf) / grid_area

                                if args.hotspots:
                                    hs = pocket.hot_spots(H)
                                    if hs is not None:
                                        hs = H >= hs
                                        hsmol2.write_scatter(D[0][hs], H[hs])
                                    else:
                                        hsmol2.write_scatter([], [])

                                volumes = []
                                for I, mol2 in zip(pocket.outer_inner(D[-1], args.io_threshold), wmol2):
                                    mol2.write_scatter(D[0][I], H[I])
                                    volumes.append(sum(I) * grid_area)
                                pockets_volume.write(('%d\t%d\t%0.1f\t%0.1f' % (window + tuple(volumes))) + os.linesep)
                            #elif args.wfull:
                            elif (wnr == 0) and (args.wfull or (not many_windows)):
                                D = pocket.distribution(paths, grid_size=grid_size, edges=edges, window=window,
                                                        pbar=pbar, map_fun=pool.imap_unordered)
                                H = (D[-1] / float(number_of_frames)) / grid_area

                                if args.hotspots:
                                    hs = pocket.hot_spots(H)
                                    mol2 = WriteMOL2(rdir + 'hotspots_full%s%s.mol2' % (ptn, args.output_suffix))
                                    if hs is not None:
                                        hs = H >= hs
                                        mol2.write_scatter(D[0][hs], H[hs])
                                    else:
                                        mol2.write_scatter([], [])
                                    del mol2

                                volumes = []
                                for I, mol2 in zip(pocket.outer_inner(D[-1], args.io_threshold),
                                                   [WriteMOL2(rdir + 'outer_full%s%s.mol2' % (ptn, args.output_suffix)),
                                                    WriteMOL2(rdir + 'inner_full%s%s.mol2' % (ptn, args.output_suffix))]):
                                    mol2.write_scatter(D[0][I], H[I])
                                    volumes.append(sum(I) * grid_area)
                                    del mol2
                                pockets_volume.write(('%d\t%d\t%0.1f\t%0.1f' % (window + tuple(volumes))) + os.linesep)
                            save_cric()
                        pool.close()
                        pool.join()

                        if wmol2 is not None:
                            for mol2 in wmol2:
                                del mol2
                        if hsmol2 is not None:
                            if args.hotspots:
                                del hsmol2
                        pockets_volume.close()

                clui.message("what it's got in its nassty little pocketses?")

            # ----------------------------------------------------------------------#
            # load mater paths data
            if args.master:
                if args.master_radius and args.master:
                    with clui.tictoc('Loading master paths data'):
                        # get stage IV options
                        options4 = config.get_stage_options(3)

                        with clui.fbm('Loading data dump from %s file' % options4.dump):
                            vda = get_vda_reader(options4.dump, mode='r')
                            result4 = vda.load()
                            mps = result4.pop(
                                'master_paths_smooth')  # FIXME: let user decide what kind of master paths are used

                        options6 = config.get_stage_options(5)

                        linmet = get_linearize_method(options6.simply_smooths)
                else:
                    # FIXME: Temporary solution for loading master paths when they are not needed due to master_radius default value.
                    mps = {}

                # ----------------------------------------------------------------------#
                # re load paths?

                if args.master_radius and args.raw_master and (not args.raw):
                    with clui.tictoc('ReLoading paths (raw data)'):
                        # get stage II options
                        options2 = config.get_stage_options(1)
                        with clui.fbm('Loading data dump from %s file' % options2.dump):
                            vda = get_vda_reader(options2.dump, mode='r')
                            result2 = vda.load()
                        with clui.fbm('Preparing raw data'):
                            paths = result2.pop('paths')
                            if args.raw_singl:
                                for p in paths:
                                    p.discard_singletons(singl=args.raw_singl)
                            paths = [p for p in paths if len(p.frames)]
                            if paths_types:
                                paths = [p for p in paths if p.name in paths_types]

                # ----------------------------------------------------------------------#
                # master paths profiles

                limit_ctypes = [ct.strip() for ct in args.master_ctypes.split(' ')]
                if limit_ctypes != [""]:
                    clui.message('Limiting master paths data to %s ctypes.' % (' '.join(limit_ctypes)))
                    for ctk in list(mps.keys()):
                        if str(ctk) not in limit_ctypes:
                            mps.pop(ctk)

                if args.master_radius and len(mps) == 0:
                    clui.message("No master paths data.")

                if args.eng_pro and len(mps) and ref is not None:
                    with clui.tictoc('Master paths profiles calculation'):

                        limit_ctypes = [ct.strip() for ct in args.master_ctypes.split(' ')]
                        if limit_ctypes != [""]:
                            for ctk in list(mps.keys()):
                                if str(ctk) not in limit_ctypes:
                                    mps.pop(ctk)

                        W = args.windows
                        WS = args.wsize

                        many_windows = W > 1 or (W == 1 and (WS is not None))
                        many_windows = many_windows and WS < Reader.number_of_frames(onelayer=True)

                        for wnr, window in enumerate(
                                pocket.windows(Reader.number_of_frames(onelayer=True), windows=W, size=WS)):

                            number_of_frames = (window[-1] - window[0])
                            if Reader.sandwich_mode:
                                number_of_frames *= Reader.number_of_layers()

                            if wnr or args.wfull:  # or not many_windows:
                                for ctype, mp in mps.items():
                                    if isinstance(mp,dict):
                                        logger.warning("Pond cannot yet handle MasterPaths calculated with separate_master option.")
                                        logger.warning("MasterPaths for %s skip." % ctype)
                                        continue

                                    fname_window = ""
                                    fname_window_single = ""
                                    fname = str(ctype).replace(':', '-')
                                    if not wnr and args.wfull:
                                        fname += 'full'
                                    elif wnr and many_windows:
                                        fname_window = '_W%d' % wnr
                                        fname_window_single = "_WX"

                                    window_name = fname if fname else wnr
                                    with clui.pbar(mess="Calculate energy for window {}".format(window_name),
                                                   maxval=int(number_of_frames)) as pbar:

                                        pool = Pool(processes=optimal_threads.threads_count)

                                        mode = 'w'
                                        if wnr > 1:
                                            mode = 'a'

                                        centers, ids = linmet(mp.coords_cont, ids=True)

                                        D = pocket.sphere_density_raw(Reader.iterate(), args.ref_mol,
                                                                      centers, args.master_radius, pool,
                                                                      window=window, pbar=pbar)

                                        H = D / float(number_of_frames) / (4. / 3. * np.pi * float(args.master_radius) ** 3)
                                        if ref:
                                            H = -k * args.temp * np.log(H) - ref

                                        with WriteMOL2(rdir + "mp_%s%s_radius%s%s.mol2" % (fname, fname_window_single, ptn, args.output_suffix),
                                                       mode=mode) as mol2:
                                            mol2.write_connected(centers, H)

                                        with open(rdir + "mp_%s%s_radius.dat" % (fname, fname_window), 'w') as dat:
                                            dat.write('len\tE' + os.linesep)
                                            L = np.hstack((0., np.cumsum(traces.diff(centers))))
                                            for l, E in zip(L, H):
                                                dat.write('%f\t%f%s' % (l, E, os.linesep))
                                    pool.close()
                                    pool.join()
                                save_cric()

            # ----------------------------------------------------------------------#
            # Paths
            # ----------------------------------------------------------------------#
            # Extracting

            def path_id_formatter(id):
                """
                Change tuple values to string with values separated by colon. Otherwise return string representation.
                :param id: tuple
                :return: string
                """
                if isinstance(id, tuple):
                    return ":".join([str(x) for x in id])
                else:
                    return str(id)

            if args.extract_path:
                if not args.raw_path:
                    options3 = config.get_stage_options(2)
                    with clui.fbm("Loading data dump from {} file.".format(options3.dump)):
                        vda = get_vda_reader(options3.dump, mode="r")
                        result3 = vda.load()

                        paths = result3.pop("spaths")
                else:
                    options2 = config.get_stage_options(1)
                    with clui.fbm('Loading raw data dump from %s file' % options2.dump):
                        vda = get_vda_reader(options2.dump, mode='r')
                        result2 = vda.load()

                        paths = result2.pop("paths")

                with clui.fbm("Extracting path {} coordinates.".format(args.extract_path)):
                    try:
                        path = next(path for path in paths if path_id_formatter(path.id) == str(args.extract_path))

                    except StopIteration:
                        path = None

                    if path:
                        if not args.raw_path:
                            soptions = result3.pop("soptions")
                            soptions = namedtuple('Options', list(soptions.keys()))(*list(soptions.values()))
                            smooth_method = get_smooth_method(soptions)

                            coords = path.get_coords_cont(smooth_method)
                        else:
                            coords = path.coords

                        with open(rdir + args.output_file, "w") as csvfile:
                            csv_writer = csv.writer(csvfile)

                            for coord in coords:
                                csv_writer.writerow(coord)
                    else:
                        clui.message("Path with ID {} does not exists.".format(args.extract_path))

                if path:
                    clui.message("Path saved to {}".format(args.output_file))

            # ----------------------------------------------------------------------#
            # Calculating energy profile for path

            if args.eng_pro and (args.path_id or args.path_file):
                with clui.tictoc('Loading paths data'):

                    if not args.raw_path:
                        options3 = config.get_stage_options(2)
                        with clui.fbm("Loading data dump from {} file.".format(options3.dump)):
                            vda = get_vda_reader(options3.dump, mode="r")
                            result3 = vda.load()

                            paths = result3.pop("spaths")
                    else:
                        options2 = config.get_stage_options(1)
                        with clui.fbm('Loading raw data dump from %s file' % options2.dump):
                            vda = get_vda_reader(options2.dump, mode='r')
                            result2 = vda.load()

                            paths = result2.pop("paths")

                    if args.path_id:
                        try:
                            path = next(path for path in paths if path_id_formatter(path.id) == args.path_id)

                            if not args.raw_path:
                                if args.path_smooth:
                                    soptions = result3.pop("soptions")
                                    soptions = namedtuple('Options', list(soptions.keys()))(*list(soptions.values()))
                                    smooth_method = get_smooth_method(soptions)
                                    coords = path.get_coords_cont(smooth_method)
                                else:
                                    coords = path.get_coords_cont()
                            else:
                                coords = path.coords
                        except StopIteration:
                            coords = []
                    elif args.path_file:
                        coords = []
                        with open(args.path_file, "r") as csvfile:
                            csv_reader = csv.reader(csvfile)

                            for row in csv_reader:
                                coords.append([float(r) for r in row])

                        coords = np.array(coords, dtype=np.float32)

                    if args.raw_path and args.raw_singl:
                        with clui.fbm('Preparing raw data'):
                            for p in paths:
                                p.discard_singletons(singl=args.raw_singl)
                        paths = [p for p in paths if len(p.frames)]

                    options6 = config.get_stage_options(5)
                    linmet = get_linearize_method(options6.simply_smooths)

                    coords = linmet(coords)

                if len(coords):
                    with clui.tictoc("Path profiles calculation", stdout=sys.stdout):

                        W = args.windows
                        WS = args.wsize

                        many_windows = W > 1 or (W == 1 and (WS is not None))
                        many_windows = many_windows and WS < Reader.number_of_frames(onelayer=True)

                        path_name = args.path_id if args.path_id else ""

                        for wnr, window in enumerate(
                                pocket.windows(Reader.number_of_frames(onelayer=True), windows=W, size=WS)):

                            number_of_frames = (window[-1] - window[0])
                            if Reader.sandwich_mode:
                                number_of_frames *= Reader.number_of_layers()

                            if wnr or args.wfull:  # or not many_windows:
                                # TODO: Handling names variables need to be refactored
                                fname_window = ""
                                fname_window_single = ""
                                fname = ""
                                if not wnr and args.wfull:
                                    window_name = "full"
                                    fname += 'full'
                                elif wnr and many_windows:
                                    fname_window = '_W%d' % wnr
                                    fname_window_single = "_WX"

                                window_name = fname if fname else wnr
                                with clui.pbar(mess="Calculate energy for window {}".format(window_name),
                                               maxval=int(number_of_frames)) as pbar:

                                    pool = Pool(processes=optimal_threads.threads_count)

                                    mode = 'w'
                                    if wnr > 1:
                                        mode = 'a'

                                    D = pocket.sphere_density_raw(Reader.iterate(), args.ref_mol,
                                                                  coords, args.path_radius, pool,
                                                                  window=window, pbar=pbar)

                                    H = D / float(number_of_frames) / (
                                            4. / 3. * np.pi * float(args.path_radius) ** 3)
                                    if ref:
                                        H = -k * args.temp * np.log(H) - ref

                                    with WriteMOL2(
                                            rdir + "path%s_%s%s_radius%s%s.mol2" % (
                                                    path_name, fname, fname_window_single, ptn, args.output_suffix),
                                            mode=mode) as mol2:
                                        mol2.write_connected(coords, H)

                                    with open(rdir + "path%s_%s%s_radius.dat" % (path_name, fname, fname_window),
                                              'w') as dat:
                                        dat.write('len\tE' + os.linesep)
                                        L = np.hstack((0., np.cumsum(traces.diff(coords))))
                                        for l, E in zip(L, H):
                                            dat.write('%f\t%f%s' % (l, E, os.linesep))

                                    pool.close()
                                    pool.join()
                            save_cric()
                else:
                    clui.message("Path with ID {} does not exists.".format(args.path_id))

            Reader.reset()
            with gzip.open(rdir + ('pond_meta%s%s.json' % (ptn, args.output_suffix)), mode='wt', compresslevel=9) as f:
                json.dump(results_meta, f)
                # TODO: consider usage of IterEncoder - move it to aquaduct/apps/data.py module

            # ----------------------------------------------------------------------#

    except BaseException:
        clui.emit_tvtb_to_file_in_root_logger(exc_info())
        raise
