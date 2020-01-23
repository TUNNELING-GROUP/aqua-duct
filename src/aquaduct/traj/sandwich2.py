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

################################################################################
# rest of imports

import re
from os.path import splitext
from os import pathsep
from collections import OrderedDict, namedtuple

import numpy as np

from MDAnalysis.topology.core import guess_atom_element

from aquaduct.utils.helpers import SmartRange, SmartRangeIncrement
from aquaduct.geom.convexhull import SciPyConvexHull, are_points_within_convexhull
from aquaduct.utils.helpers import arrayify, create_tmpfile, tupleify
from aquaduct.utils.maths import make_default_array
from aquaduct.apps.data import GCS, CRIC
from aquaduct.utils.maths import defaults
from aquaduct.utils.helpers import ion

################################################################################
# import or create memory decorator

if GCS.cachedir:
    from joblib import Memory
    memory_cache = Memory(cachedir=GCS.cachedir,
                          verbose=0)
    # mmap have to be switched off, otherwise smoothing does not work properly
    # memory_cache = Memory(cachedir=GCS.cachedir, mmap_mode='r', verbose=0)
    memory = memory_cache.cache
elif GCS.cachemem:
    from aquaduct.utils.helpers import memory_in_memory as memory
else:
    from aquaduct.utils.helpers import noaction as memory


################################################################################
# Trajectory data

class MetaReader(object):

    baguette = 'baguette'
    sandwich = 'sandwich'
    shortbread = 'shortbread'

    def __init__(self, topology=[], trajectory=[],
                 mode='baguette',
                 window=slice(None),
                 threads=1,
                 engine='mda'):

        self.topology = topology
        self.trajectory = trajectory
        self.mode = mode
        self.window = window
        self.threads = threads
        self.engine = engine

    @property
    def number_of_readers(self):
        if self.mode == self.sandwich:
            return len(self.topology)
        return self.threads

    def iter_readers(self):
        if self.mode == self.sandwich:
            for top in self.topology:



################################################################################

class ReaderTraj(object):

    def __init__(self, topology=None, trajectory=[],
                 window=slice(None)):
        self.topology = topology
        self.trajectory = trajectory
        self.window = window

        self.trajectory_object = self.open()

    def open(self):
        # should return any object that can be further used to parse trajectory
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def close(self):
        # should close trajectory reader in self.trajectory_object
        # WARNING: This method has to be carefully implemented because it is used by __del__ and
        #          should not emmit any error messages. This is subjet of change.
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def __del__(self):
        return self.close()


################################################################################
# raw opener

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

def open_raw(topology, trajectory, engine):
    if engine == 'mda':
        return open_raw_mda(topology,trajectory)

################################################################################

class ReaderTrajViaMDA(ReaderTraj):

    def open(self):
        return open_raw_mda(self.topology,self.trajectory)

    def close(self):
        if hasattr(self, 'trajectory_object'):
            if hasattr(self.trajectory_object, 'trajectory'):
                if hasattr(self.trajectory_object.trajectory, 'close'):
                    self.trajectory_object.trajectory.close()


################################################################################

# flag sandwich as imported
GCS.sandwich_import = True
