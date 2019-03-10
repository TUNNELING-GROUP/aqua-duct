# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2016-2018  Tomasz Magdziarz, Alicja Płuciennik, Michał Stolarczyk <info@aquaduct.pl>
# Copyright (C) 2019  Tomasz Magdziarz <info@aquaduct.pl>
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

import numpy as np

# from scipy.spatial import ConvexHull as SciPyConvexHull
from scipy.spatial import ConvexHull

from aquaduct.utils.helpers import uniqify
from aquaduct.geom.traces import vector_change_len


################################################################################
# SciPy ConvexHull improvements

class SciPyConvexHull(ConvexHull):

    def __init__(self, points, inflate=None, *args, **kwargs):
        super(SciPyConvexHull, self).__init__(points, *args, **kwargs)
        # at this stage we have chull
        if inflate:
            inflate = float(inflate)
            center = np.mean(self.vertices_points, 0)
            # redo chull but inflate its vertices
            new_points = (vector_change_len(p - center, inflate) + center for p in points)
            super(SciPyConvexHull, self).__init__(list(new_points), *args, **kwargs)

    @property
    def vertices_ids(self):
        return self.vertices.tolist()

    @property
    def vertices_points(self):
        return self.points[self.vertices]

    @property
    def facets(self):
        for simp in self.simplices:
            yield np.vstack((self.points[simp], self.points[simp][0]))

    @property
    @uniqify
    def edges(self):
        for simp in self.simplices:
            n = len(simp)
            for nr, point in enumerate(simp):
                if nr < n - 1:
                    yield point, simp[nr + 1]
                else:
                    yield point, simp[0]

    @property
    def simplices_vertices(self):
        for simp in self.simplices:
            yield [self.vertices_ids.index(s) for s in simp]

    def point_within(self, point):
        # checks if point is within convexhull
        # build new convex hull with this point and vertices
        new_hull = SciPyConvexHull(np.vstack((point, self.vertices_points)))
        # return 0 not in _unfold_vertices_list(new_hull.vertices)
        return 0 not in new_hull.vertices_ids


################################################################################

# primitives to check if point(s) is within chull

def is_point_within_convexhull(point_chull):
    # This is helper function to check if point is within convex hull.
    # point_chull is a tuple where first element holds a point and the second is a ConvexHull object.
    return point_chull[-1].point_within(point_chull[0])


# many points - new solution but soon will be deprecated

def are_points_within_convexhull(points, chull):
    vertices_points = chull.points[chull.vertices]

    promise = ((SciPyConvexHull(np.vstack((point, vertices_points))).vertices[0] != 0) for point in points)

    return np.fromiter(promise, bool)
