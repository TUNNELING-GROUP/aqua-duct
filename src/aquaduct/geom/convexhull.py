# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2016-2017  Tomasz Magdziarz, Alicja Płuciennik, Michał Stolarczyk <info@aquaduct.pl>
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

from scipy.spatial import ConvexHull as SciPyConvexHull

from aquaduct.utils.helpers import uniqify
from scipy.spatial.distance import cdist, pdist


def _vertices_ids(convexhull):
    return convexhull.vertices.tolist()


def _vertices_points(convexhull):
    # returns points of all vertices
    return convexhull.points[convexhull.vertices]


def _point_within_convexhull(convexhull, point):
    # checks if point is within convexhull
    # build new convex hull with this point and vertices
    new_hull = SciPyConvexHull(np.vstack((point, convexhull.vertices_points)))
    # return 0 not in _unfold_vertices_list(new_hull.vertices)
    return 0 not in new_hull.vertices_ids


def _facets(convexhull):
    for simp in convexhull.simplices:
        yield np.vstack((convexhull.points[simp], convexhull.points[simp][0]))


@uniqify
def _edges(convexhull):
    for simp in convexhull.simplices:
        n = len(simp)
        for nr, point in enumerate(simp):
            if nr < n - 1:
                yield point, simp[nr + 1]
            else:
                yield point, simp[0]


SciPyConvexHull.vertices_ids = property(_vertices_ids)
SciPyConvexHull.vertices_points = property(_vertices_points)
SciPyConvexHull.point_within = _point_within_convexhull
SciPyConvexHull.facets = property(_facets)
SciPyConvexHull.edges = property(_edges)


def is_point_within_convexhull(point_chull):
    # This is helper function to check if point is within convex hull.
    # point_chull is a tuple where first element holds a point and the second is a ConvexHull object.
    return point_chull[-1].point_within(point_chull[0])


'''
class ConvexHull(object):
    def __init__(self,points):
        self.points = points

        self.chulls = self.collection_of_chulls(self.points)

    def get_normal_chull(self,points):
        return SciPyConvexHull(points)

    def find_vertices_connected_by_longest_edge(self,chull):
        lengths = [float(pdist(chull.points[edge,:])) for edge in chull.edges]
        return chull.edges[np.argmax(lengths)]

    def split_chull(self,chull):
        edges = self.find_vertices_connected_by_longest_edge(chull)
        vertices = chull.vertices_ids
        vertices.pop(vertices.index[edges[0]])
        yield self.get_normal_chull(chull.points[vertices])
        vertices = chull.vertices_ids
        vertices.pop(vertices.index[edges[1]])
        yield self.get_normal_chull(chull.points[vertices])

    def collection_of_chulls(self,points):
        return list(self.split_chull(self.get_normal_chull(points)))
'''
