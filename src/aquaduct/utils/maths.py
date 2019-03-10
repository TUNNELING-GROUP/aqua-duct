# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2016-2018  Tomasz Magdziarz, Alicja Płuciennik, Michał Stolarczyk <info@aquaduct.pl>
# Copyright (C) 2019  Tomasz Magdziarz <info@aquaduct.pl>
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
    """
    Default types that are enforced in :class:`numpy.ndarray` objects.

    .. note::

        It is used only througt :attr:`defaults` instance.

    """
    float_default = np.float64
    int_default = np.int64
    int_type = np.int8


defaults = NumpyDefaultsStorageTypes()
"""
Instance of :class:`~NumpyDefaultsStorageTypes` to store default values.
"""


def make_default_array(array_like):
    """
    :param array_like: Array like object
    :return: Array with dtype set to :attr:`NumpyDefaultsStorageTypes.float_default`.
    """
    if isinstance(array_like, np.ndarray):
        return array_like.astype(defaults.float_default)
    return np.array(array_like).astype(defaults.float_default)


class MemMap(namedtuple('MemMap', 'filename dtype shape')):
    """
    Provides simple convenience wrapper for :func:`numpy.memmap`.
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
    Convenience class for handling :class:`numpy.ndarray` and :class:`numpy.core.memmap.memmap` objects in a transparent way.
    """

    def __init__(self, filename=None, dtype=None, shape=None):
        """
        :param filename str: Optional name of the file to store memory mapped object.
        :param dtype: Optional dtype of array, if `None` default value of :class:`NumpyDefaultsStorageTypes.float_default` is used.
        :param shape: Shape of the array.

        If no `filename` is given then regular :class:`numpy.ndarray` is created with :func:`numpy.zeros`.
        Otherwise :class:`~MemMap` object is created.
        """
        assert shape is not None, "Cannot make array without shape."
        if dtype is None:
            dtype = defaults.float_default
        if filename is None:
            self.array = np.zeros(shape, dtype=dtype)
        else:
            self.array = MemMap(filename, dtype, shape)

    @property
    def isndarray(self):
        """
        :return: `True` if underlaying object is of :class:`numpy.ndarray` type.
        :rtype: bool
        """
        return isinstance(self.array)

    def readwrite(self):
        """
        :return: Array with read-write access.
        """
        if self.isndarray:
            return self.array
        return self.array.readwrite()

    def readonly(self):
        """
        :return: Array with read only access, if possible
        """
        if self.isndarray:
            return self.array
        return self.array.readonly()

    def __call__(self):
        """
        By default this calls :func:`~ArrayOrArray.readwrite`.

        :return: Array with read-write access.
        """
        return self.readwrite()
