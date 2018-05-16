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

from aquaduct.utils.helpers import uniqify
from aquaduct.geom import Sphere, do_cut_thyself

from scipy.spatial.distance import cdist, pdist, squareform


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

SciPyConvexHull.vertices_ids = property(_vertices_ids)
SciPyConvexHull.vertices_points = property(_vertices_points)
SciPyConvexHull.point_within = _point_within_convexhull
SciPyConvexHull.facets = property(_facets)
SciPyConvexHull.edges = property(_edges)
SciPyConvexHull.simplices_vertices = property(_simplices_vertices)


def is_point_within_convexhull(point_chull):
    # This is helper function to check if point is within convex hull.
    # point_chull is a tuple where first element holds a point and the second is a ConvexHull object.
    return point_chull[-1].point_within(point_chull[0])

def are_points_within_convexhull(points,chull,map_fun=None):
    are = np.zeros(len(points),dtype=bool)
    are[list(ids_points_within_convexhull_alt(points,chull,map_fun=map_fun))] = True
    return are.tolist()

def ids_points_within_convexhull_alt(points, chull, ids=None, map_fun=None):
    if map_fun is None:
        map_fun = map
    if ids is None:
        ids = np.ones(len(points)) == 1
    tree = KDTree(points)
    if len(points)>3:
        search_p = SciPyConvexHull(points).vertices_points[0]
    else:
        search_p = np.zeros(3)
    dver = squareform(pdist(chull.vertices_points)) # inter vertices distances
    output = []
    while ids.any():
        i = tree.query(search_p,k=1)[-1]
        if not ids[i]:
            i = int(np.argwhere(ids)[0])
        p = np.hstack((points[i],1))
        dmin = np.abs(np.dot(chull.equations,np.array([p]).T)) # distances to faces
        if is_point_within_convexhull((p[:3],chull)):
            dmin = dmin.min()
            output.append(i)
            ids[i] = False
            tt = tree.query_ball_point(p[:3], dmin)
            output.extend(tt)
            ids[tt] = False
        else:
            dfacets = (cdist(f[:3,:],[p[:3]]) for f,d in zip(chull.facets,dmin)) # distance p to facets points, by facets
            dfac = []
            for dfacet,simp,dm in zip(dfacets,chull.simplices_vertices,dmin):
                dmin_vertices = False
                for nr in range(3):
                    t = map(float,[dfacet[nr],dfacet[(nr+1)%3],dver[simp[nr],simp[(nr+1)%3]]])
                    if np.arccos((t[1]**2 + t[2]**2 - t[0]**2)/(2*t[1]*t[2])) > np.pi/2:
                        dmin_vertices = True
                        break
                    if np.arccos((t[0]**2 + t[2]**2 - t[1]**2)/(2*t[0]*t[2])) > np.pi/2:
                        dmin_vertices = True
                        break
                if dmin_vertices:
                    dfac.append(min(dfacet))
                else:
                    dfac.append(dm)
            dmin = min(dfac)
            tt = tree.query_ball_point(p[:3], dmin)
            ids[tt] = False
            if sum(ids)>3:
                search_p = SciPyConvexHull(points[ids, :]).vertices_points[0]
            else:
                search_p = np.zeros(3)

    return output


def ids_points_within_convexhull(points, chull, ids=None, map_fun=None):
    if map_fun is None:
        map_fun = map
    return map_fun(is_point_within_convexhull, ((oc,chull) for oc in points))

    if ids is None:
        ids = np.array(range(len(points)))
    # 0) can we calculate chull?
    if len(points) < 4:
        are_points = np.array(list(map_fun(is_point_within_convexhull, ((oc, chull) for oc in points))))
        if are_points.any():
            return ids[are_points].tolist()
        return []
    # 1) make chull
    points_chull = SciPyConvexHull(points)
    are_points_chull = np.array(list(map_fun(is_point_within_convexhull, ((oc,chull) for oc in points_chull.vertices_points))))
    # special case if all are in
    if are_points_chull.all():
        return ids.tolist()
    # 2) take all
    all_points = points_chull.vertices_points
    all_points = np.hstack((all_points,np.ones((all_points.shape[0],1)))) # add column of 1
    # 3) for each such point calculate minimal distance to faces of chull
    dmin = (np.abs(np.dot(chull.equations,all_points.T))/np.sqrt(np.matrix(all_points[:,:3]**2).sum(1)).T.A).min(0)
    # 4) for matching points construct spheres, cut thyself, sort
    within_ids = []
    if sum(are_points_chull) > sum(~are_points_chull):
        spheres = [Sphere(center,radius,nr) for nr,(center,radius) in enumerate(zip(all_points[are_points_chull,:3],dmin[are_points_chull]))]
        nrspheres = do_cut_thyself(spheres)[0] # nonredundant
        nrspheres = sorted(nrspheres,key=lambda s: -s.radius)
        # 5) find all points that are included in these spheres - they are also in chull
        d = cdist([sphe.center for sphe in nrspheres],points)
        d -= np.repeat([[sphe.radius for sphe in nrspheres]],len(points),0).T
        d = (d <= 0).any(0)
        within_ids.extend(ids[d].tolist())
        ids = ids[~d]
        points = points[~d]
        '''
        for sphe in nrspheres: # TODO: this loop can be replaced by array arithmetics
            d = cdist([sphe.center],points).flatten()
            within_ids.extend(ids[d <= sphe.radius].tolist())
            ids = ids[~(d <= sphe.radius)]
            points = points[~(d <= sphe.radius)]
        '''
        if (~are_points_chull).any():
            # 6) for non matching points onstruct spheres, cut thyself, sort
            spheres = [Sphere(center,radius,nr) for nr,(center,radius) in enumerate(zip(all_points[~are_points_chull,:3],dmin[~are_points_chull]))]
            nrspheres = do_cut_thyself(spheres)[0] # nonredundant
            nrspheres = sorted(nrspheres,key=lambda s: -s.radius)
            # 5) find all points that are included in these spheres - they are also outside chull
            d = cdist([sphe.center for sphe in nrspheres],points)
            d -= np.repeat([[sphe.radius for sphe in nrspheres]],len(points),0).T
            d = (d <= 0).any(0)
            ids = ids[~d]
            points = points[~d]
            '''
            for sphe in nrspheres: # TODO: this loop can be replaced by array arithmetics
                d = cdist([sphe.center],points).flatten()
                ids = ids[~(d <= sphe.radius)]
                points = points[~(d <= sphe.radius)]
            '''
    else:
        if (~are_points_chull).any():
            # 6) for non matching points onstruct spheres, cut thyself, sort
            spheres = [Sphere(center,radius,nr) for nr,(center,radius) in enumerate(zip(all_points[~are_points_chull,:3],dmin[~are_points_chull]))]
            nrspheres = do_cut_thyself(spheres)[0] # nonredundant
            nrspheres = sorted(nrspheres,key=lambda s: -s.radius)
            # 5) find all points that are included in these spheres - they are also outside chull
            d = cdist([sphe.center for sphe in nrspheres],points)
            d -= np.repeat([[sphe.radius for sphe in nrspheres]],len(points),0).T
            d = (d <= 0).any(0)
            ids = ids[~d]
            points = points[~d]
            '''
            for sphe in nrspheres: # TODO: this loop can be replaced by array arithmetics
                d = cdist([sphe.center],points).flatten()
                ids = ids[~(d <= sphe.radius)]
                points = points[~(d <= sphe.radius)]
            '''
        if (are_points_chull).any():
            spheres = [Sphere(center,radius,nr) for nr,(center,radius) in enumerate(zip(all_points[are_points_chull,:3],dmin[are_points_chull]))]
            nrspheres = do_cut_thyself(spheres)[0] # nonredundant
            nrspheres = sorted(nrspheres,key=lambda s: -s.radius)
            # 5) find all points that are included in these spheres - they are also in chull
            d = cdist([sphe.center for sphe in nrspheres],points)
            d -= np.repeat([[sphe.radius for sphe in nrspheres]],len(points),0).T
            d = (d <= 0).any(0)
            within_ids.extend(ids[d].tolist())
            ids = ids[~d]
            points = points[~d]
            '''
            for sphe in nrspheres: # TODO: this loop can be replaced by array arithmetics
                d = cdist([sphe.center],points).flatten()
                within_ids.extend(ids[d <= sphe.radius].tolist())
                ids = ids[~(d <= sphe.radius)]
                points = points[~(d <= sphe.radius)]
            '''
    if len(ids):
            return within_ids + ids_points_within_convexhull(points, chull, ids, map_fun=map_fun)
    return within_ids

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
        B = np.random.randn(10000, 3)*10
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
