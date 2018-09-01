# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2018  Tomasz Magdziarz <info@aquaduct.pl>
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

from aquaduct.traj.paths import GenericPaths

import numpy as np
from itertools import izip,imap

from scipy import spatial

def get_spc(sp,window=None):
    '''
    :param sp: Single path like object or Generic path.
    :param tuple window: Optional frames window.
    :rtype: numpy.ndarray
    :return: Coordinates of path; to be used in pocket calculation.
    '''
    if isinstance(sp,GenericPaths):
        if window is None:
            return sp.coords
        i = (np.array(sp.frames)>=window[0]) & (np.array(sp.frames)<=window[1])
        return get_spc(sp)[i]
    else:
        if window is None:
            return sp.coords_cont
        i = (np.array(sp.paths_cont)>=window[0]) & (np.array(sp.paths_cont)<=window[1])
        return get_spc(sp)[i]


def find_minmax(spaths,pbar=None):
    '''
    :param list spaths: List of single like path objects.
    :param pbar: Optional progress object providing next() method.
    :rtype: 2 element tuple of numpy.ndarray each of shape (3,)
    :return: Minimal and maximal boundaries of coordinates used in pocket calulations of spaths.
    '''
    minc = np.array([float('inf')]*3)
    maxc = np.array([float('-inf')]*3)
    for sp in spaths:
        sp_minc = get_spc(sp).min(0)
        minc_i = minc > sp_minc
        minc[minc_i] = sp_minc[minc_i]
        sp_maxc = get_spc(sp).max(0)
        maxc_i = maxc < sp_maxc
        maxc[maxc_i] = sp_maxc[maxc_i]
        if pbar:
            pbar.next()
    minc = np.floor(minc)
    maxc = np.ceil(maxc)
    return minc,maxc

def find_minmax_single(sp):
    coords = get_spc(sp)
    return coords.min(0),coords.max(0)

def find_minmax_map(spaths,pbar=None,map_fun=None):
    '''
    :param list spaths: List of single like path objects.
    :param pbar: Optional progress object providing next() method.
    :rtype: 2 element tuple of numpy.ndarray each of shape (3,)
    :return: Minimal and maximal boundaries of coordinates used in pocket calulations of spaths.
    '''
    minc = np.array([float('inf')]*3)
    maxc = np.array([float('-inf')]*3)
    if map_fun is None:
        map_fun = imap
    for sp_minc,sp_maxc in map_fun(find_minmax_single,spaths):
        minc_i = minc > sp_minc
        minc[minc_i] = sp_minc[minc_i]
        maxc_i = maxc < sp_maxc
        maxc[maxc_i] = sp_maxc[maxc_i]
        if pbar:
            pbar.next()
    minc = np.floor(minc)
    maxc = np.ceil(maxc)
    return minc,maxc

def find_edges(spaths,grid_size=1.,pbar=None,map_fun=None):
    '''
    :param list spaths: List of single like path objects.
    :param float grid_size: Size of grid cell in A.
    :param pbar: Optional progress object providing next() method.
    :rtype: list of numpy.ndarrays
    :return: Edges of bins of grid spanning all submited paths.
    '''
    return [np.linspace(mi,ma,int((ma-mi)*grid_size)+1) for mi,ma in zip(*find_minmax_map(spaths,pbar=pbar,map_fun=map_fun))]

class distribution_worker(object):
    def __init__(self,edges=None,window=None):
        self.edges = edges
        self.window = window
    def __call__(self,sp):
        return np.histogramdd(get_spc(sp,window=self.window),bins=self.edges)[0]

def distribution(spaths,grid_size=1.,edges=None,window=None,pbar=None,map_fun=None):
    '''
    :param list spaths: List of single like path objects.
    :param float grid_size: Size of grid cell in A.
    :param list of numpy.ndarrays edges: Edges of bins of grid spanning all submited paths.
    :param tuple window: Optional frames window.
    :param pbar: Optional progress object providing next() method.
    :rtype tuple of numpy.ndarrays
    :return: Coordinates of pocket and number of points.
    '''
    maxc = np.array(map(max,edges))
    minc = np.array(map(min,edges))
    H = np.zeros(map(int, (maxc - minc) * grid_size))
    if map_fun is None:
        map_fun = map
    map_worker = distribution_worker(edges=edges,window=window)
    for h in map_fun(map_worker,spaths):
        H += h
        if pbar:
            pbar.next()
    mg = [ee[:-1]+(1./(grid_size+1)) for ee in edges]
    x,y,z = np.meshgrid(*mg,indexing='ij')
    pocket = H > 0
    return np.vstack((x[pocket],y[pocket],z[pocket])).T,H[pocket]

def outer_inner(H):
    '''
    :param numpy.ndarray H: Pocket distribution.
    :return: Indices of outer and inner pocket.
    :rtype: tuple of numpy.ndarray
    '''
    OI = H / H[H>0].mean()
    return OI<1,OI>=1

def windows(frames,windows=None,size=None):
    yield (0.,frames - 1.) # full window
    if windows:
        if size:
            begs = np.linspace(0, frames- size, windows)
            ends = np.array([b+size-1 for b in begs])
            if ends[-1] > frames - 1:
                ends[-1] = frames - 1
        else:
            begs = np.linspace(0, frames, windows + 1)[:-1]
            ends = np.linspace(-1,frames-1,windows+1)[1:]
            if ends[-1] < frames - 1:
                ends[-1] = frames - 1
        for b,e in zip(begs,ends):
            yield np.floor(b),np.floor(e)


def sphere_radii(spaths,centers=None,radii=None,window=None,pbar=None,map_fun=None):

    H = np.zeros(len(centers),dtype=np.int32)

    if map_fun is None:
        map_fun = imap

    for sp in spaths:
        coords = get_spc(sp,window=window)
        D = spatial.distance.cdist(coords,centers)
        for nr,radius in enumerate(radii):
            H[nr] += np.count_nonzero(D[:,nr]<=radius)
        if pbar:
            pbar.next()
    return H

class sphere_radius_worker(object):
    def __init__(self,window,centers,radius):
        self.window = window
        self.centers = centers
        self.radius = radius
    def __call__(self,sp):
        coords = get_spc(sp,window=self.window)
        D = spatial.distance.cdist(coords,self.centers) <= self.radius
        return np.count_nonzero(D,0)

def sphere_radius(spaths,centers=None,radius=2.,window=None,pbar=None,map_fun=None):

    H = np.zeros(len(centers),dtype=np.int32)

    if map_fun is None:
        map_fun = imap
    map_worker = sphere_radius_worker(window,centers,radius)
    for h in map_fun(map_worker,spaths):
        H += h
        if pbar:
            pbar.next()
    return H

def hot_spots(H):
    return hot_spots_his(H)

def hot_spots_his(H,bins=(10,21)):
    bn = []
    for b in xrange(*bins):
        his = np.histogram(H,bins=b)
        if 0 in his[0]:
            i = int(np.argwhere(his[0]==0)[0])
            bn.append((b,sum(his[0][i:])))
    if len(bn):
        i = np.argmax(np.array(bn)[:,1])
        b = bn[i][0]
        his = np.histogram(H,bins=b)
        hs = np.argwhere(his[0]==0)
        return his[1][hs][0]
