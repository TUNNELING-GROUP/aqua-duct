import numpy as np

from pyhull.convex_hull import ConvexHull

def _unfold_vertices_list(vertices):
    vv = []
    [vv.extend(v) for v in vertices]
    vv = list(set(vv))
    vv.sort()
    return vv

def _vertices_ids(convexhull):
    return _unfold_vertices_list(convexhull.vertices)

def _vertices_points(convexhull):
    # returns points of all vertices
    return np.array([convexhull.points[v] for v in convexhull.vertices_ids])

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
        yield np.vstack((simp.coords,simp.coords[0]))

ConvexHull.facets = property(_facets)
