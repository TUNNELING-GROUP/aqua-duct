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

from scipy.spatial import ConvexHull as SciPyConvexHull

from aquaduct.utils.helpers import uniqify
from aquaduct.geom import Sphere, do_cut_thyself

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


def are_points_within_convexhull(points,chull,ids=None,map_fun=None):
    if map_fun is None:
        map_fun = map
    if ids is None:
        ids = np.array(range(len(points)))
    # 1) make chull
    points_chull = SciPyConvexHull(points)
    are_points_chull = np.array(list(map_fun(is_point_within_convexhull, ((oc,chull) for oc in points_chull.vertices_points))))
    # special case if all are in
    if are_points_chull.all():
        return range(len(points))
    # 2) take matching only
    are_points = points_chull.vertices_points[are_points_chull]
    are_points = np.hstack((are_points,np.ones((are_points.shape[0],1))))
    # 3) for each such point calculate minimal distance to faces of chull
    dmin = (np.abs(np.dot(chull.equations,are_points.T))/np.sqrt(np.matrix(are_points[:,:3]**2).sum(1)).T.A).min(0)
    # 4) starting from the biggest what is cluded in spheres is also within
    spheres = [Sphere(center,radius,nr) for nr,(center,radius) in enumerate(zip(are_points[:,:3],dmin))]
    nrspheres = do_cut_thyself(spheres)[0]
    nrspheres = sorted(nrspheres,key=lambda s: -s.radius)
    within_ids = []
    for sphe in nrspheres:
        d = cdist([sphe.center],points).flatten()
        within_ids.extend(ids[d <= sphe.radius].tolist())
        ids = ids[~(d <= sphe.radius)]
        points = points[~(d <= sphe.radius)]
    if len(within_ids) and len(points):
        return within_ids+are_points_within_convexhull(points,chull,ids,map_fun=map_fun)

    return within_ids

if __name__ == "__main__":
    import numpy as np

    A = np.random.randn(100, 3)
    B = np.random.randn(100, 3)
    chull = SciPyConvexHull(A)

    old_way = np.argwhere(map(lambda p: is_point_within_convexhull((p,chull)),B)).flatten()

    new_way = sorted(are_points_within_convexhull(B,chull))
