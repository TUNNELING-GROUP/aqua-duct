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

# from matplotlib.colors import colorConverter

import scipy.spatial.distance as distance

from aquaduct.traj.paths import GenericPathTypeCodes as gptc
from aquaduct.traj.paths import PathTypesCodes as ptc
from aquaduct.visual.cmaps import default as default_cmap
from aquaduct.utils.helpers import zip_zip, is_number, lind

_cl2rgba = {'b': (0.0, 0.0, 1.0, 1.0),
            'c': (0.0, 0.75, 0.75, 1.0),
            'g': (0.0, 0.5, 0.0, 1.0),
            'k': (0.0, 0.0, 0.0, 1.0),
            'm': (0.75, 0, 0.75, 1.0),
            'r': (1.0, 0.0, 0.0, 1.0),
            'y': (0.75, 0.75, 0, 1.0)}


def euclidean(A,B):
    return distance.cdist(A, B, 'euclidean')

def cityblock(A,B):
    return distance.cdist(A, B, 'cityblock')

# cc = lambda c, alpha=1.0: colorConverter.to_rgb(c)

def cc_safe(c):
    # color converter
    if c in 'rgbcmyk':
        c = _cl2rgba[c]
    c = tuple(c)
    assert len(c) in [3, 4], 'Color should be given either as one letter (rgbcmyk) or as rgb or rgba vector.'
    for e in c:
        assert is_number(e), 'Color vector has to be specified as numbers.'
        assert 0 <= e <= 1, 'Color vector elements have to be in range of 0 to 1.'
    return c[:3]


def cc(c):
    # color converter faster
    if c in 'rgbcmyk':
        c = _cl2rgba[c]
    return c[:3]


_dcc_is = ptc.path_in_code + gptc.scope_name
_dcc_cc = ptc.path_object_code + gptc.object_name
_dcc_cs = ptc.path_object_code + gptc.scope_name
_dcc_os = ptc.path_out_code + gptc.scope_name

_dcc_i = ptc.path_in_code
_dcc_c = ptc.path_object_code
_dcc_o = ptc.path_out_code

_default_color_codes = {_dcc_is: 'r',
                        _dcc_cc: 'g',
                        _dcc_cs: 'y',
                        _dcc_os: 'b',
                        _dcc_i: 'r',
                        _dcc_c: 'g',
                        _dcc_o: 'b'}

default_color_codes = _default_color_codes


def color_codes(code, custom_codes=None):
    if custom_codes is None:
        return default_color_codes[code]
    else:
        return custom_codes[code]


def get_cmap(size):
    return [e[0][0] for e in zip_zip(default_cmap, N=size)]


class ColorMapDistMap(object):

    grey = (0.5, 0.5, 0.5, 1)

    def __init__(self):
        self.cmap = default_cmap
        self.cmap = self.__do_cadex()

    def distance(self,E1,E2):
        # E1 and E2 are RGBs matrices
        D = []
        for e1 in E1:
            DD = []
            for e2 in E2:
                DD.append(self.color_distance(e1,e2))
            D.append(DD)
        return np.array(D)

    @staticmethod
    def color_distance(e1,e2):
        # e1 and e2 are colors defs of rgb components in this order
        # Taken from http://www.compuphase.com/cmetric.htm
        # scale them to 0-255 range
        def to0255(E):
            return [int(e*255) for e in E]
        e1 = to0255(e1)
        e2 = to0255(e2)

        rmean = (e1[0]+e2[0])/2
        r = e1[0]-e2[0]
        g = e1[1]-e2[1]
        b = e1[2]-e2[2]
        return np.sqrt((((512 + rmean) * r * r) >> 8) + 4 * g * g + (((767 - rmean) * b * b) >> 8))

    def __do_cadex(self):
        m = len(self.cmap) # number of objects
        k = len(self.cmap)
        # indices
        mi = [] # model
        ti = range(m) # test
        # FIRST OBJECT
        # get distance to mean object
        Xm = np.array(self.cmap).mean(axis=0) # mean object
        Xm.shape = (1,Xm.shape[0]) # fix shape
        Dm = self.distance(Xm,np.array(self.cmap))
        D = self.distance(np.array(self.cmap),np.array(self.cmap))
        # min value will be the first object
        #mi.append(ti.pop(Dm.argmin()))
        mi.append(ti.pop(-1))
        if k > 1:
            # SECOND OBJECT
            Xm = np.array(self.cmap)[mi[-1]]
            Xm.shape = (1,Xm.shape[0]) # fix shape
            Dm = self.distance(Xm,np.array(self.cmap)[ti,:])
            # max value will be the first object
            mi.append(ti.pop(Dm.argmax()))
            if k > 2:
                while (len(mi)<k):
                    Dm = D[:,ti][mi,:]
                    mi.append(ti.pop(Dm.min(axis=0).argmax()))
        return lind(self.cmap,mi+ti)

    def __call__(self, node):
        if 0 < node:
            # cycle...
            return self.cmap[(node-1) % len(self.cmap)][:3]
            #return self.cmap[int(np.round(self.cm_size * f_like(node)))][:3]
        # return grey otherwise
        return self.grey[:3]


def f_like(n):
    if n == 1:
        return 0.0
    if n == 2:
        return 0.5
    n -= 1
    order = np.floor(np.log(n) / np.log(2))
    parts = 2 ** order
    current = n - parts
    return 0.5 / parts + 1. / parts * current
