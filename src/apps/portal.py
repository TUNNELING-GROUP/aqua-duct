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
This is a magic portal detector.
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

if __name__ == "__main__":

    from aquaduct.utils import clui

    with clui.tictoc('Magic Portal calculations'):

        ############################################################################
        # argument parsing
        import argparse
        from aquaduct import version_nice as aquaduct_version_nice

        description_version = '''Aquaduct library version %s''' % (aquaduct_version_nice(),)
        description = '''Magic Portal (debugging tool, not for regular use)'''

        parser = argparse.ArgumentParser(description=description,
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)

        parser.add_argument("-t", action="store", dest="top", required=False, default=None,
                            help="Topology file.")
        parser.add_argument("-p", action="store", dest="tra", required=False, default=None,
                            help="Trajectory file(s).")
        parser.add_argument("-a", action="store", dest="areas", required=False, default=None,
                            help="Portal areas, semicolon separated.")

        args = parser.parse_args()

        areas = args.areas.split(';')

        from aquaduct.traj.sandwich import Reader, Window
        from os import pathsep

        Reader(args.top, [trj.strip() for trj in args.tra.split(pathsep)], Window(None, None, 1))  # trajectory reader

        from aquaduct.geom.convexhull import SciPyConvexHull

        print("frame", end=' ')
        for nr, area in enumerate(areas):
            print(("a%da" % (nr + 1)), ("a%dv" % (nr + 1)), end=' ')
        print("")

        for traj_reader in Reader.iterate():
            traj_reader = traj_reader.open()
            for frame in traj_reader.iterate_over_frames():
                print(frame, end=' ')
                for area in areas:
                    chull = traj_reader.parse_selection(area).chull()
                    print(chull.area, chull.volume, end=' ')
                print("")
