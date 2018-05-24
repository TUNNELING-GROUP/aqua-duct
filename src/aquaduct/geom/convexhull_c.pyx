# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2016-2018  Tomasz Magdziarz, Alicja Płuciennik, Michał Stolarczyk <info@aquaduct.pl>
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

import numpy as np
cimport cython
cimport numpy as np

from scipy.spatial import ConvexHull

ctypedef np.float_t DTYPE

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cpdef np.ndarray are_points_within_convexhull(np.ndarray[DTYPE, ndim=2] points, chull):
    cdef np.ndarray[DTYPE, ndim=2] vertices_points = chull.points[chull.vertices]
    cdef np.ndarray[DTYPE, ndim=1] point

    promise = (np.vstack((point,vertices_points)) for point in points)
    #promise = ((ConvexHull(np.vstack((point,vertices_points))).vertices[0] != 0) for point in points)

    return np.fromiter(promise,bool)



