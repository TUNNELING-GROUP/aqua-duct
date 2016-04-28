'''
Created on Dec 4, 2015

@author: tljm
'''

from aqueduct.geom.convexhull_scipy import ConvexHull

def is_point_within_convexhull(point,chull):
    return chull.point_within(point)
