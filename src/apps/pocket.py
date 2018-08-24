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

from aquaduct.apps.data import GCS, load_cric


if __name__ == "__main__":

    from aquaduct.utils import clui

    with clui.tictoc('What have I got in my pocket?'):

        ############################################################################
        # argument parsing
        import argparse
        from aquaduct import version_nice as aquaduct_version_nice

        description_version = '''Aquaduct library version %s''' % (aquaduct_version_nice(),)
        description = '''What have I got in my pocket?'''

        parser = argparse.ArgumentParser(description=description,
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)

        parser.add_argument("-c", action="store", dest="config_file", required=False, help="Config file filename.")
        parser.add_argument("--max-frame", action="store", dest="max_frame", type=int, required=False,
                            help="Maximal number of frame.")
        parser.add_argument("--min-frame", action="store", dest="min_frame", type=int, required=False,
                            help="Minimal number of frame.")
        parser.add_argument("--step-frame", action="store", dest="step_frame", type=int, required=False,
                            help="Frames step.")
        parser.add_argument("--sandwich", action="store_true", dest="sandwich", required=False,
                            help="Sandwich mode for multiple trajectories.")
        parser.add_argument("--cache-dir", action="store", dest="cachedir", type=str, required=False,
                            help="Directory for coordinates caching.")
        parser.add_argument("--window-full", action="store_true", dest="wfull", required=False,
                            help="Return full windows.")
        parser.add_argument("--windows", action="store", dest="windows", type=int, required=True,
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

        args = parser.parse_args()

        ############################################################################

        clui.message('What have I got in my pocket?')

        ############################################################################
        # cache dir & netcdf
        GCS.cachedir = args.cachedir
        load_cric()

        from aquaduct.traj.sandwich import Reader, Window
        from aquaduct.apps.valvecore import ValveConfig, valve_load_config
        from aquaduct.apps.valvedata import get_vda_reader

        ############################################################################
        # config
        config = ValveConfig()  # config template
        valve_load_config(args.config_file, config)

        # get global options
        goptions = config.get_global_options()

        # Maximal frame checks
        frames_window = Window(args.min_frame, args.max_frame, args.step_frame)

        # Open trajectory reader
        # trajectories are split with os.pathsep
        from os import pathsep

        Reader(goptions.top, [trj.strip() for trj in goptions.trj.split(pathsep)],
               window=frames_window,
               sandwich=args.sandwich)  # trajectory reader

        ############################################################################
        # load spaths

        # get stage III options
        options3 = config.get_stage_options(2)

        with clui.fbm('Loading data dump from %s file' % options3.dump):
            vda = get_vda_reader(options3.dump, mode='r')
            result3 = vda.load()
            spaths = result3.pop('spaths')

        ############################################################################
        # run pockets?

        import numpy as np

        ref = None
        # caclulte n_z for reference
        if args.ref:
            for traj_reader in Reader.iterate():
                traj_reader = traj_reader.open()
                for frame in traj_reader.iterate():
                    ref_sel = traj_reader.parse_selection(args.ref)
                    com = ref_sel.center_of_mass()
                    break
                break
            del traj_reader

            with clui.pbar(Reader.number_of_frames(onelayer=True),mess="Calculating density in reference:") as pbar:
                ref = 0.
                for traj_reader in Reader.iterate():
                    traj_reader = traj_reader.open()
                    ref_sel = '(%s) and (point %f %f %f %f)' % ((args.ref_mol,)+tuple(map(float,com))+(args.ref_radius,))
                    for frame in traj_reader.iterate():
                        ref += len(list(traj_reader.parse_selection(ref_sel).residues().ids()))
                        pbar.next()

            ref /= Reader.number_of_frames(onelayer=True)
            ref /= 4./3.*np.pi*(float(args.ref_radius)**3)

            # B constant
            k = 0.0019872041 # kcal/mol/K
            ref = -k*args.temp*np.log(ref)


        from aquaduct.geom import pocket
        from aquaduct.traj.dumps import WriteMOL2

        W = args.windows
        WS = args.wsize
        grid_size = 1.
        grid_area = grid_size ** 3
        with clui.pbar(len(spaths) * (1 + W + 1), mess='Calculating pockets:') as pbar:
            edges = pocket.find_edges(spaths, grid_size=grid_size, pbar=pbar)
            number_of_frames = Reader.number_of_frames(onelayer=True)
            if WS is None:
                WSf = float(number_of_frames/float(W))
            else:
                WSf = float(WS)
            wmol2 = [WriteMOL2('outer.mol2'), WriteMOL2('inner.mol2')]
            for wnr, window in enumerate(pocket.windows(number_of_frames, windows=W, size=WS)):
                if wnr:
                    D = pocket.distribution(spaths, grid_size=grid_size, edges=edges, window=window, pbar=pbar)
                    H = (D[-1] / WSf)/grid_area
                    if ref:
                        H = -k*args.temp*np.log(H) - ref
                    for I, mol2 in zip(pocket.outer_inner(D[-1]), wmol2):
                        mol2.write_scatter(D[0][I], H[I])
                else:
                    D = pocket.distribution(spaths, grid_size=grid_size, edges=edges, window=window, pbar=pbar)
                    H = (D[-1] / float(number_of_frames))/grid_area
                    if ref:
                        H = -k*args.temp*np.log(H) - ref
                    for I, mol2 in zip(pocket.outer_inner(D[-1]), [WriteMOL2('outer_full.mol2'), WriteMOL2('inner_full.mol2')]):
                        mol2.write_scatter(D[0][I], H[I])
                        del mol2

            for mol2 in wmol2:
                del mol2

        ############################################################################

        clui.message("what it's got in its nassty little pocketses?")

        ############################################################################
