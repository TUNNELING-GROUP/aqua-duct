# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2018-2019  Tomasz Magdziarz, Michał Banas <info@aquaduct.pl>
# Copyright (C) 2020  Tomasz Magdziarz <info@aquaduct.pl>
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


from aquaduct import logger
from aquaduct.utils.helpers import version_parser
from os.path import splitext
import re

import MDAnalysis as mda

################################################################################
# Check MDAnalysis version

def mda_ver():
    return version_parser(mda.__version__)


# FIXME: do it according to user options
if mda_ver() < version_parser('0.16'):
    logger.error('Unsupported MDAnalysis version %s; should be 0.16.2 or > 0.19.', mda.__version__)
    raise NotImplementedError('Unsupported MDAnalysis version %s; should be 0.16.2 or > 0.19.' % mda.__version__)

if mda_ver() < version_parser('0.17') and mda_ver() < version_parser('0.20'):
    logger.warning('Unsupported MDAnalysis version %s.', mda.__version__)
    logger.warning('Trying to mitigate potential problems by setting `use_periodic_selections = False`.')
    mda.core.flags['use_periodic_selections'] = False


ENGINE_MDA = 'mda'

available_engines = [ENGINE_MDA]


################################################################################
# mda

mda_available_formats = {re.compile('(nc|NC)'): 'nc',
                         re.compile('(prmtop|parmtop|top|PRMTOP|PARMTOP|TOP)'): 'PRMTOP',
                         re.compile('(dcd|DCD)'): 'LAMMPS',
                         re.compile('(psf|PSF)'): 'psf',
                         re.compile('(pdb|PDB)'): 'pdb',
                         re.compile('(crd|CRD)'): 'crd',
                         re.compile('(xtc|XTC)'): 'XTC'}

def open_raw_mda(topology, trajectory):
    topology_ext = splitext(topology)[1][1:]
    for afk in list(mda_available_formats.keys()):
        if afk.match(topology_ext):
            topology_ext = mda_available_formats[afk]
            break
    trajectory_ext = splitext(trajectory[0])[1][1:]
    for afk in list(mda_available_formats.keys()):
        if afk.match(trajectory_ext):
            trajectory_ext = mda_available_formats[afk]
            break
    return mda.Universe(topology,
                        trajectory,
                        topology_format=topology_ext,
                        format=trajectory_ext)

################################################################################

def open_raw(topology, trajectory, engine):
    if engine == ENGINE_MDA:
        return open_raw_mda(topology,trajectory)

################################################################################
