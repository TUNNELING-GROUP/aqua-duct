#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2018  Tomasz Magdziarz
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



################################################################################

class count_ref(object):

    def __init__(self,pbar,ref_sel):

        self.pbar = pbar # queue
        self.ref_sel = ref_sel

    def __call__(self,traj_reader):
        traj_reader = traj_reader.open()
        ref = 0.
        for frame in traj_reader.iterate():
            ref += len(list(traj_reader.parse_selection(self.ref_sel).residues().ids()))
            self.pbar.put(1) # report one frame a time
        return ref


################################################################################


if __name__ == "__main__":

    from aquaduct.utils import clui

    with clui.tictoc('What have I got in my pocket?'):

        #----------------------------------------------------------------------#
        # argument parsing

        import argparse
        from aquaduct import version_nice as aquaduct_version_nice

        description_version = '''Aquaduct library version %s''' % (aquaduct_version_nice(),)
        description = '''What have I got in my pocket?'''

        parser = argparse.ArgumentParser(description=description,
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)

        parser.add_argument("-c", action="store", dest="config_file", required=False, help="Config file filename.")
        parser.add_argument("-t", action="store", dest="threads", required=False, default=None,
                            help="Limit Aqua-Duct calculations to given number of threads.")
        parser.add_argument("-r", action="store", dest="results_dir", required=False, help="Path to results directory",default="",type=str)
        parser.add_argument("--max-frame", action="store", dest="max_frame", type=int, required=False,
                            help="Maximal number of frame.")
        parser.add_argument("--min-frame", action="store", dest="min_frame", type=int, required=False,
                            help="Minimal number of frame.")
        parser.add_argument("--step-frame", action="store", dest="step_frame", type=int, required=False,
                            help="Frames step.")
        parser.add_argument("--raw", action="store_true", dest="raw", required=False,
                            help="Use raw data from paths instead of single paths.")
        parser.add_argument("--raw-master", action="store_true", dest="raw_master", required=False,
                            help="Use raw data from paths instead of single paths, only in master paths calculations.")
        parser.add_argument("--raw-discard-singletons", action="store", dest="raw_singl", required=False, type=int, default=1,
                            help="Discard short scope only segments from raw data.")
        parser.add_argument("--sandwich", action="store_true", dest="sandwich", required=False,
                            help="Sandwich mode for multiple trajectories.")
        parser.add_argument("--cache-dir", action="store", dest="cachedir", type=str, required=False,
                            help="Directory for coordinates caching.")
        parser.add_argument("--window-full", action="store_true", dest="wfull", required=False,
                            help="Return full window if windows is used.")
        parser.add_argument("--windows", action="store", dest="windows", type=int, required=False, default=1,
                            help="Number of windows to calculate.")
        parser.add_argument("--wsize", action="store", dest="wsize", type=int, required=False, default=None,
                            help="Size of window in frames.")
        parser.add_argument("--reference", action="store", dest="ref", type=str, required=False, default=None,
                            help="Selection of reference in the first frame of trajectory.")
        parser.add_argument("--reference-radius", action="store", dest="ref_radius", type=float, required=False, default=2.,
                            help="Radius of reference.")
        parser.add_argument("--reference-mol", action="store", dest="ref_mol", type=str, required=False, default='resname WAT',
                            help="Selection of reference molecules.")
        parser.add_argument("--temperature", action="store", dest="temp", type=float, required=False, default=300.,
                            help="Simulation temperature.")
        parser.add_argument("--gsize", action="store", dest="grid_size", type=float, required=False, default=1.,
                            help="Size of grid's cells.")
        parser.add_argument("--pockets", action="store_true", dest="pockets", required=False,
                            help="Calculate pockets.")
        parser.add_argument("--hotspots", action="store_true", dest="hotspots", required=False,
                            help="Calculates hotspots if pockets are calculated.")
        parser.add_argument("--master-radius", action="store", dest="master_radius", type=float, required=False,
                            help="Calculate profiles for master paths with giwen radius.")
        parser.add_argument("--master-ctypes", action="store", dest="master_ctypes", type=str, required=False, default="",
                            help="Limit calculations to given ctypes.")
        #parser.add_argument("--master-radii", action="store_true", dest="master_radii", required=False,
        #                    help="Calculate profiles for master paths using width as radii.")

        args = parser.parse_args()


        #----------------------------------------------------------------------#
        # cache dir
        from aquaduct.apps.data import GCS, load_cric, save_cric

        GCS.cachedir = args.cachedir
        load_cric()

        from aquaduct.traj.sandwich import Reader, Window
        from aquaduct.apps.valve.core import valve_load_config, ValveConfig
        from aquaduct.apps.valve.data import get_vda_reader
        from aquaduct.utils.multip import optimal_threads
        from os import pathsep
        import numpy as np
        from multiprocessing import Pool,Manager
        from aquaduct.geom import pocket
        from aquaduct.traj.dumps import WriteMOL2
        from aquaduct.apps.valve.helpers import get_linearize_method
        from aquaduct.geom import traces
        from itertools import izip
        import os
        from aquaduct.geom.smooth import SavgolSmooth


        #----------------------------------------------------------------------#
        # results dir
        rdir = ""
        if len(args.results_dir):
            if not os.path.exists(args.results_dir):
                os.makedirs(args.results_dir)
            rdir = args.results_dir
        rdir += os.path.sep

        #----------------------------------------------------------------------#
        # config
        config = ValveConfig()  # config template
        valve_load_config(args.config_file, config)

        # get global options
        goptions = config.get_global_options()


        if args.threads is None:
            optimal_threads.threads_count = optimal_threads.cpu_count + 1
        else:
            optimal_threads.threads_count = int(args.threads)
        clui.message("Number of threads Pond is allowed to use: %d" % optimal_threads.threads_count)
        if (1 < optimal_threads.threads_count < 3) or (optimal_threads.threads_count - 1 > optimal_threads.cpu_count):
            clui.message(
                "Number of threads is not optimal; CPU count reported by system: %d" % optimal_threads.cpu_count)
        # because it is used by mp.Pool it should be -1???
        if optimal_threads.threads_count > 1:
            optimal_threads.threads_count -= 1
            clui.message("Main process would use 1 thread.")
            clui.message("Concurent calculations would use %d threads." % optimal_threads.threads_count)


        # Maximal frame checks
        frames_window = Window(args.min_frame, args.max_frame, args.step_frame)

        # Open trajectory reader
        # trajectories are split with os.pathsep

        Reader(goptions.top, [trj.strip() for trj in goptions.trj.split(pathsep)],
               window=frames_window,
               sandwich=args.sandwich,
               threads=optimal_threads.threads_count) # trajectory reader

        #----------------------------------------------------------------------#
        # load paths

        if (args.pockets or args.master_radius) and not (args.master_radius and not args.raw and args.raw_master):
            if not args.raw:
                # get stage III options
                options3 = config.get_stage_options(2)
                with clui.fbm('Loading data dump from %s file' % options3.dump):
                    vda = get_vda_reader(options3.dump, mode='r')
                    result3 = vda.load()
                paths = result3.pop('spaths')
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

        #----------------------------------------------------------------------#
        # reference value

        ref = None
        # calculate n_z for reference
        if args.ref:
            for traj_reader in Reader.iterate():
                traj_reader = traj_reader.open()
                for frame in traj_reader.iterate():
                    ref_sel = traj_reader.parse_selection(args.ref)
                    com = ref_sel.center_of_mass()
                    break
                break
            Reader.reset()

            goal = Reader.number_of_frames(onelayer=False)
            with clui.pbar(goal,mess="Calculating density in the reference area:") as pbar:
                ref = []
                ref_sel = '(%s) and (point %f %f %f %f)' % ((args.ref_mol,)+tuple(map(float,com))+(args.ref_radius,))
                #ref_sel = '(%s) and (sphzone %f %s)' % (args.ref_mol, args.ref_radius, args.ref)
                manager = Manager()
                pbar_queue = manager.Queue()
                pool = Pool(processes=optimal_threads.threads_count)

                r = pool.map_async(
                    count_ref(pbar_queue,ref_sel),
                    Reader.iterate(),
                    callback=ref.extend)

                progress = 0
                for p in iter(pbar_queue.get,None):
                    progress += p
                    pbar.next(step=p)
                    if progress == goal:
                        break
                r.wait()
                pool.close()
                pool.join()

            ref = float(sum(ref))

            clui.message('Reference number of molecules: %d [molecules].' % int(ref))

            ref /= float(Reader.number_of_frames(onelayer=False)) # if sandwich mean value will be correct

            ref /= 4./3.*np.pi*(float(args.ref_radius)**3)
            clui.message('Reference density: %0.4f [molecules/A^3].' % ref)

            # Boltzmann k constant
            #k = 0.0019872041 # kcal/mol/K
            k = 0.0083144621 # kJ/mol/K

            ref = -k*args.temp*np.log(ref)
            clui.message('Reference correction: %0.4f [kJ/mol/K].' % ref)

        #----------------------------------------------------------------------#
        # calculate pockets

        if args.pockets:

            clui.message('What have I got in my pocket?')

            W = args.windows
            WS = args.wsize
            grid_size = args.grid_size
            grid_area = grid_size ** 3
            with clui.pbar(len(paths) * (1 + W + int(args.wfull)), mess='Calculating pockets:') as pbar:
                pool = Pool(processes=optimal_threads.threads_count)
                edges = pocket.find_edges(paths, grid_size=grid_size, pbar=pbar,map_fun=pool.imap_unordered)
                number_of_frames = Reader.number_of_frames(onelayer=True)
                if WS is None:
                    WSf = float(number_of_frames/float(W))
                else:
                    WSf = float(WS)
                wmol2 = [WriteMOL2(rdir+'outer.mol2'), WriteMOL2(rdir+'inner.mol2')]
                if args.hotspots:
                    hsmol2 = WriteMOL2(rdir+'hotspots.mol2')
                for wnr, window in enumerate(pocket.windows(Reader.number_of_frames(onelayer=True), windows=W, size=WS)):
                    number_of_frames = (window[-1]-window[0])
                    if Reader.sandwich_mode:
                        number_of_frames *= Reader.number_of_layers()

                    if wnr:
                        D = pocket.distribution(paths, grid_size=grid_size, edges=edges, window=window, pbar=pbar, map_fun=pool.imap_unordered)
                        H = (D[-1] / WSf)/grid_area
                        if args.hotspots:
                            hs = pocket.hot_spots(H)
                            if hs is not None:
                                hs = H >= hs
                        if ref:
                            H = -k*args.temp*np.log(H) - ref
                        if args.hotspots:
                            if hs is not None:
                                hsmol2.write_scatter(D[0][hs], H[hs])
                            else:
                                hsmol2.write_scatter([], [])
                        for I, mol2 in zip(pocket.outer_inner(D[-1]), wmol2):
                            mol2.write_scatter(D[0][I], H[I])
                    elif args.wfull:
                        D = pocket.distribution(paths, grid_size=grid_size, edges=edges, window=window, pbar=pbar, map_fun=pool.imap_unordered)
                        H = (D[-1] / float(number_of_frames))/grid_area
                        if ref:
                            H = -k*args.temp*np.log(H) - ref
                        for I, mol2 in zip(pocket.outer_inner(D[-1]), [WriteMOL2(rdir+'outer_full.mol2'), WriteMOL2(rdir+'inner_full.mol2')]):
                            mol2.write_scatter(D[0][I], H[I])
                            del mol2
                    save_cric()
                pool.close()
                pool.join()

                for mol2 in wmol2:
                    del mol2
                if args.hotspots:
                    del hsmol2


            clui.message("what it's got in its nassty little pocketses?")

        #----------------------------------------------------------------------#
        # load mater paths data

        if args.master_radius:

            # get stage IV options
            options4 = config.get_stage_options(3)

            with clui.fbm('Loading data dump from %s file' % options4.dump):
                vda = get_vda_reader(options4.dump, mode='r')
                result4 = vda.load()
                mps = result4.pop('master_paths_smooth')

            options6 = config.get_stage_options(5)

            linmet = get_linearize_method(options6.simply_smooths)

        #----------------------------------------------------------------------#
        # re load paths?

        if args.master_radius and args.raw_master and (not args.raw):
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

        #----------------------------------------------------------------------#
        # master paths profiles

        limit_ctypes = [ct.strip() for ct in args.master_ctypes.split(' ')]
        if limit_ctypes != [""]:
            clui.message('Limiting master paths data to %s ctypes.' % (' '.join(limit_ctypes)))
            for ctk in mps.keys():
                if str(ctk) not in limit_ctypes:
                    mps.pop(ctk)

        if args.master_radius and len(mps) == 0:
            clui.message("No master paths data.")

        if args.master_radius and len(mps):

            limit_ctypes = [ct.strip() for ct in args.master_ctypes.split(' ')]
            if limit_ctypes != [""]:
                for ctk in mps.keys():
                    if str(ctk) not in limit_ctypes:
                        mps.pop(ctk)


            W = args.windows
            WS = args.wsize

            many_windows = W > 1 or (W==1 and  (WS is not None))
            many_windows = many_windows and WS < Reader.number_of_frames(onelayer=True)

            pbar_len = len(paths)*len(mps) * (W + int(args.wfull))

            pool = Pool(processes=optimal_threads.threads_count)

            with clui.pbar(pbar_len, mess='Calculating master paths profiles:') as pbar:

                #number_of_frames = Reader.number_of_frames(onelayer=False)
                #window = pocket.windows(Reader.number_of_frames(onelayer=True), windows=W, size=WS).next() # only full window

                for wnr,window in enumerate(pocket.windows(Reader.number_of_frames(onelayer=True), windows=W, size=WS)):

                    #print many_windows,wnr,window

                    number_of_frames = (window[-1]-window[0])
                    if Reader.sandwich_mode:
                        number_of_frames *= Reader.number_of_layers()

                    if wnr or args.wfull: # or not many_windows:

                        for ctype,mp in mps.iteritems():

                            fname_window = ""
                            fname_window_single = ""
                            fname = str(ctype).replace(':','-')
                            if not wnr and args.wfull:
                                fname += '_full'
                            elif wnr and many_windows:
                                fname_window = '_W%d' % wnr
                                fname_window_single = "_WX"

                            mode = 'w'
                            if wnr > 1:
                                mode = 'a'

                            centers, ids = linmet(mp.coords_cont, ids=True)

                            #print "fname: '%s' fname_window: '%s' mode: '%s'" % (fname,fname_window,mode)


                            D = pocket.sphere_radius(paths,centers=centers,radius=args.master_radius,window=window,pbar=pbar,map_fun=pool.imap_unordered)
                            #print number_of_frames
                            H = D / float(number_of_frames) / (4./3. * np.pi * float(args.master_radius)**3)
                            if ref:
                                H = -k*args.temp*np.log(H) - ref
                            with WriteMOL2(rdir+"mp_%s%s_radius.mol2" % (fname,fname_window_single),mode=mode) as mol2:
                                mol2.write_connected(centers,H)

                            with open(rdir+"mp_%s%s_radius.dat" % (fname,fname_window),'w') as dat:
                                dat.write('len\tE'+os.linesep)
                                L = np.hstack((0.,np.cumsum(traces.diff(centers))))
                                for l,E in izip(L,H):
                                    dat.write('%f\t%f%s' % (l,E,os.linesep))
                        save_cric()

            pool.close()
            pool.join()

        #----------------------------------------------------------------------#

        Reader.reset()

        #----------------------------------------------------------------------#
