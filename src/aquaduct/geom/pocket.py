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


import numpy as np

def get_spc(sp,window=None):
    '''
    :param MacromolPath sp: Single path like object.
    :param tuple window: Optional frames window.
    :rtype: numpy.ndarray
    :return: Coordinates of path; to be used in pocket calculation.
    '''
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

def find_edges(spaths,grid_size=1.,pbar=None):
    '''
    :param list spaths: List of single like path objects.
    :param float grid_size: Size of grid cell in A.
    :param pbar: Optional progress object providing next() method.
    :rtype: list of numpy.ndarrays
    :return: Edges of bins of grid spanning all submited paths.
    '''
    return [np.linspace(mi,ma,int((ma-mi)*grid_size)+1) for mi,ma in zip(*find_minmax(spaths,pbar=pbar))]

def distribution(spaths,grid_size=1.,edges=None,window=None,pbar=None):
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
    for sp in spaths:
        H += np.histogramdd(get_spc(sp,window=window),bins=edges)[0]
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



'''

W = 5
with clui.pbar(len(spaths)*(1+W+1),mess='Calculating pockets:') as pbar:

    number_of_frames = Reader.number_of_frames(onelayer=True)


    A = 1./1 # grid size in A
    AS = 1 / A

    minc = np.array([float('inf')]*3)
    maxc = np.array([float('-inf')]*3)
    for sp in spaths:
        sp_minc = get_spc(sp).min(0)
        minc_i = minc > sp_minc
        minc[minc_i] = sp_minc[minc_i]
        sp_maxc = get_spc(sp).max(0)
        maxc_i = maxc < sp_maxc
        maxc[maxc_i] = sp_maxc[maxc_i]
        pbar.next()
    minc = np.floor(minc)
    maxc = np.ceil(maxc)
    e = [np.linspace(mi,ma,int((ma-mi)*AS)+1)  for mi,ma in zip(minc,maxc)]



    for wnr,window in enumerate([(0,number_of_frames-1)] + zip(np.linspace(0,number_of_frames-1,W+1)[:-1],np.linspace(0,number_of_frames-1,W+1)[1:])):

        H = np.zeros(map(int, (maxc - minc) * AS))
        for sp in spaths:
            H += np.histogramdd(get_spc(sp,window=tuple(window)),bins=e)[0]
            pbar.next()

        mg = [ee[:-1]+(1./(AS+1)) for ee in e]
        x,y,z = np.meshgrid(*mg,indexing='ij')
        pocket = H > H[H>0].mean()
        pocket = H > 0
        H /= float(number_of_frames)
        HH = H / H[H>0].mean()
        H /= H[H>0].mean()


        for fiona in np.linspace(0,1,6)[:-1]:
            pocket_fiona = (HH > fiona) & (HH < fiona + 1./5)
            pocket = pocket_fiona
            pdb_name = 'grid_window%d_core%0.1f.pdb' % (wnr,fiona)
            with WritePDB(pdb_name,scale_bf=1.) as wpdb:
                wpdb.write_scatter(np.vstack((x[pocket],y[pocket],z[pocket])).T,H[pocket])
        fiona = 1.
        pocket_fiona = HH > fiona
        pocket = pocket_fiona
        pdb_name = 'grid_window%d_core%0.1f.pdb' % (wnr,fiona)
        with WritePDB(pdb_name,scale_bf=1.) as wpdb:
            wpdb.write_scatter(np.vstack((x[pocket],y[pocket],z[pocket])).T,H[pocket])

'''
