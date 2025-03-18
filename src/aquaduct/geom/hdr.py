# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2018-2019  Tomasz Magdziarz <info@aquaduct.pl>
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

'''
Module performs HDR 2D calculations only with Gaussian Kerneld Density Estimator
as impelemented in :mod:`scipy.stats`.
'''

import numpy as np
from scipy.stats import gaussian_kde
from aquaduct.geom.pca import PCA, Polarize
from aquaduct.geom.traces import square_area


class HDR(object):
    def __init__(self, X, points=10, expand_by=1., center_of_system=None):
        # expand_by in A is added to grid boudaries to expand it a bit

        self.center_of_system = center_of_system

        # calcualte PCA
        self.pca = PCA(preprocess=Polarize(center=center_of_system, rvar=0.01))
        self.pca.build(X)

        # get kernel for PC1 and PC2
        self.kernel = gaussian_kde(self.pca.T[:, [0, 1]].T)

        # find x and y ranges
        # xyr = np.array([self.pca.T[:,:2].max(0),self.pca.T[:,:2].min(0)]).mean(0) # center of range
        xyr = (self.pca.T[:, :2].max(0) - self.pca.T[:, :2].min(0) + expand_by) / 2.  # 1/2 expanded range
        self.xyr = np.array([np.array([self.pca.T[:, :2].max(0), self.pca.T[:, :2].min(0)]).mean(0) - xyr,
                             np.array([self.pca.T[:, :2].max(0), self.pca.T[:, :2].min(0)]).mean(0) + xyr])

        # desired number of points in the grid for each of two dimensions per A

        # total number of cells is points^2
        xy = self.xyr[1] - self.xyr[0]
        points = xy.mean() * points

        b = (((points ** 2) * xy[1]) / (xy[1] * (xy[0] ** 2))) ** 0.5
        a = xy[0] / xy[1] * b

        self.points = (xy * np.array([a, b]))

        # calcualte spatial span
        self.X, self.Y = np.mgrid[self.pca.T.min(0)[0]:self.pca.T.max(0)[0]:(self.points[0] * 1j),
                         self.pca.T.min(0)[1]:self.pca.T.max(0)[1]:(self.points[1] * 1j)]
        self.positions = np.vstack([self.X.ravel(), self.Y.ravel()])

        self.values = self.kernel(self.positions)

        self.Z = np.reshape(self.values, self.X.shape)

    @property
    def cell_area(self):
        return self.cell_dimensions.prod()

    @property
    def cell_dimensions(self):
        return (self.pca.T[:, :2].max(0) - self.pca.T[:, :2].min(0)) / self.points

    def _fraction2indices(self, fraction=0.9):
        # indices are returned in the order of values, i.e. after ravel.
        return self.values > self.Z.max() * (1 - fraction)

    def area_naive(self, fraction=0.9):
        '''
        This is naive implementation. It calculates are in the space of prepocessed variables.
        Transformation to normal Cartesian space may not be trivial as details mey depends on the
        actual location  of points.

        :param float fraction: Fraction of cluster for which are should be estimated.
        :rtype: float
        :return: Area of cluster calculated in the preprocessed space.
        '''
        return self.cell_area * (self._fraction2indices(fraction=fraction)).sum()

    def _square_me(self, position):
        return self.pca.undo(position - self.cell_dimensions / 2 * np.array([[1, 1], [1, -1], [-1, -1], [-1, 1]]),
                             pc=[0, 1])

    def area(self, fraction=0.9):
        '''
        :param float fraction: Fraction of cluster for which are should be estimated.
        :rtype: float
        :return: Area of cluster calculated in the Cartesian space.
        '''
        i = self._fraction2indices(fraction=fraction)
        return sum((square_area(*tuple(sq)) for sq in map(self._square_me, self.positions.T[i])))
