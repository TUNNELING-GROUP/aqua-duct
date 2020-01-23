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
                pass


################################################################################


################################################################################

# flag sandwich as imported
GCS.sandwich_import = True

