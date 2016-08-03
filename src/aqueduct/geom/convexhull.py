# -*- coding: utf-8 -*-

import numpy as np

from scipy.spatial import ConvexHull

def _vertices_ids(convexhull):
    return convexhull.vertices.tolist()

def _vertices_points(convexhull):
    # returns points of all vertices
    return convexhull.points[convexhull.vertices]

ConvexHull.vertices_ids = property(_vertices_ids)
ConvexHull.vertices_points = property(_vertices_points)

def _point_within_convexhull(convexhull,point):
    # checks if point is within convexhull
    # build new convex hull with this point and vertices
    new_hull = ConvexHull(np.vstack((point,convexhull.vertices_points)))
    #return 0 not in _unfold_vertices_list(new_hull.vertices)
    return 0 not in new_hull.vertices_ids

ConvexHull.point_within = _point_within_convexhull

def _facets(convexhull):
    for simp in convexhull.simplices:
        yield np.vstack((convexhull.points[simp],convexhull.points[simp][0]))

ConvexHull.facets = property(_facets)

def is_point_within_convexhull(point_chull):
    return point_chull[-1].point_within(point_chull[0])
