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
from scipy.spatial import KDTree
from scipy.spatial.distance import cdist, pdist, squareform

from aquaduct.utils.helpers import uniqify
from aquaduct.utils.multip import optimal_threads

################################################################################
# SciPy ConvexHull improvements

# helper functions


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


def _simplices_vertices(convexhull):
    for simp in convexhull.simplices:
        yield [convexhull.vertices_ids.index(s) for s in simp]

# add them to the class

SciPyConvexHull.vertices_ids = property(_vertices_ids)
SciPyConvexHull.vertices_points = property(_vertices_points)
SciPyConvexHull.point_within = _point_within_convexhull
SciPyConvexHull.facets = property(_facets)
SciPyConvexHull.edges = property(_edges)
SciPyConvexHull.simplices_vertices = property(_simplices_vertices)

################################################################################
# primitives to check if point(s) is within chull

def is_point_within_convexhull(point_chull):
    # This is helper function to check if point is within convex hull.
    # point_chull is a tuple where first element holds a point and the second is a ConvexHull object.
    return point_chull[-1].point_within(point_chull[0])

# many points - new solution but soon will be deprecated

def are_points_within_convexhull(points,chull,map_fun=None,sane_huge=10000):
    points = np.array(list(points))
    if map_fun is None:
        map_fun = map
        n = len(points)
    else:
        n = max(1, len(points) / optimal_threads.threads_count)
    if len(points) > sane_huge:
        return np.hstack(map_fun(are_points_within_convexhull_core,
                                 ((points[i:i + n], chull) for i in xrange(0, len(points), n))))
    return np.array(map(lambda p: is_point_within_convexhull((p,chull)),points))

def are_points_within_convexhull_core(points_chull):
    points,chull = points_chull
    ids = np.ones(len(points)) == 1
    output = np.ones(len(points)) == 0
    tree = KDTree(points)
    if len(points)>3:
        search_p = SciPyConvexHull(points).vertices_points[0]
    else:
        search_p = np.zeros(3)
    dver = squareform(pdist(chull.vertices_points)) # inter vertices distances
    search_p_change = True
    N = len(ids)/10.
    while sum(ids)>N:
        i = tree.query(search_p,k=1)[-1]
        if not ids[i]:
            i = int(np.argwhere(ids)[0])
        p = np.hstack((points[i],1))
        if is_point_within_convexhull((p[:3],chull)):
            dmin = np.abs(np.dot(chull.equations, np.array([p]).T))  # distances to faces
            dmin = dmin.min()
            output[i] = True
            ids[i] = False
            tt = tree.query_ball_point(p[:3], dmin)
            output[tt]=True
            ids[tt] = False
            if search_p_change:
                search_p_change = False
                search_p = chull.vertices_points.mean(0)
        else:
            dfacets = (cdist(f[:3,:],[p[:3]]) for f in chull.facets) # distance p to facets points, by facets
            dfac = []
            for dfacet,simp,equa in zip(dfacets,chull.simplices_vertices,chull.equations):
                dmin_vertices = 0
                for nr in range(3):
                    t = map(float,[dfacet[nr],dfacet[(nr+1)%3],dver[simp[nr],simp[(nr+1)%3]]])
                    if np.arccos((t[1]**2 + t[2]**2 - t[0]**2)/(2*t[1]*t[2])) > np.pi/2:
                        dmin_vertices += 1
                    elif np.arccos((t[0]**2 + t[2]**2 - t[1]**2)/(2*t[0]*t[2])) > np.pi/2:
                        dmin_vertices += 1
                    if dmin_vertices > 1: break
                if dmin_vertices > 1:
                    dfac.append(min(dfacet))
                else:
                    _dm = np.abs((equa*p).sum())
                    dfac.append(_dm)
            dmin = min(dfac)
            tt = tree.query_ball_point(p[:3], dmin)
            ids[tt] = False
            search_p_change = True
            if sum(ids)>3:
                search_p = SciPyConvexHull(points[ids, :]).vertices_points[0]
            else:
                search_p = np.zeros(3)
    if ids.any(): # fall back to is_point...
        output[[int(i) for i in np.argwhere(ids) if is_point_within_convexhull((points[i,:],chull))]]=True

    return output

################################################################################

if __name__ == "__main__":
    import numpy as np


    def plot_sphere(point,radius,ax,color='b',prec=30):
        x = np.cos(np.linspace(0,2*np.pi,prec))*radius
        y = np.sin(np.linspace(0,2*np.pi,prec))*radius
        ax.plot(x+point[0],y+point[1],point[2],zdir='z',color=color)
        ax.plot(x+point[0],y+point[2],point[1],zdir='y',color=color)
        ax.plot(x+point[1],y+point[2],point[0],zdir='x',color=color)


    def plot_chull(chull,ax,color='g'):
        [ax.plot(f[:,0],f[:,1],f[:,2],color=color) for f in chull.facets]




    for nr in xrange(100):

        A = np.random.randn(100, 3)
        B = np.random.randn(10000, 3)
        chull = SciPyConvexHull(A)

        old_way = np.array(map(lambda p: is_point_within_convexhull((p,chull)),B))
        new_way = np.array(are_points_within_convexhull(B, chull))
        print nr
        if (old_way != new_way).any():
            print "discrepancy!"
            dis = np.argwhere(~(old_way==new_way)).flatten()
            print dis
            for d in dis:
                print d,'old',old_way[d],'new',new_way[d]
            if len(dis) > 1:
                break
