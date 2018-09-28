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
from collections import namedtuple

class NumpyDefaultsStorageTypes(object):
    float_default = np.float64
    int_default = np.int64
    int_type = np.int8

defaults = NumpyDefaultsStorageTypes()


def make_default_array(array_like):
    if isinstance(array_like, np.ndarray):
        return array_like.astype(defaults.float_default)
    return np.array(array_like).astype(defaults.float_default)


class MemMap(namedtuple('MemMap','filename dtype shape')):
    """
    Provides simple convenience wrapper for :python:`numpy.memmap`.
    """
    def readonly(self):
        """
        :return:  Memory map object in 'r' mode.
        :rtype: numpy.core.memmap.memmap
        """
        return np.memmap(self.filename, mode='r', dtype=self.dtype, shape=self.shape)
    def readwrite(self):
        """
        :return:  Memory map object in 'r+' mode.
        :rtype: numpy.core.memmap.memmap
        """
        return np.memmap(self.filename, mode='r+', dtype=self.dtype, shape=self.shape)

class ArrayOrArray(object):
    """
    Convenience class for handling :python:`numpy.ndarray` and :python:`numpy.core.memmap.memmap` objects in a transparent way.
    """

    def __init__(self,filename=None,dtype=None,shape=None):
        assert shape is not None, "Cannot make array without shape."
        if dtype is None:
            dtype = defaults.float_default
        if filename is None:
            self.array = np.zeros(shape,dtype=dtype)
        else:
            self.array = MemMap(filename,dtype,shape)

    @property
    def isndarray(self):
        return isinstance(self.array)

    def readwrite(self):
        if self.isndarray:
            return self.array
        return self.array.readwrite()

    def readonly(self):
        if self.isndarray:
            return self.array
        return self.array.readonly()

    def __call__(self):
        return self.readwrite()

#class ArrayOrArray(namedtuple('ArrayOrArray',))