'''
Created on Dec 4, 2015

@author: tljm
'''

from aqueduct.geom.convexhull_scipy import ConvexHull

def is_point_within_convexhull(point_chull):
    return point_chull[-1].point_within(point_chull[0])
