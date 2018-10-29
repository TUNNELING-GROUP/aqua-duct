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

'''
Module performs HDR 2D calculations only with Gaussian Kerneld Density Estimator
as impelemented in :mod:`scipy.stats`.
'''

import numpy as np
from scipy.stats import gaussian_kde
from aquaduct.geom.pca import PCA

from matplotlib import _contour

class HDR(object):
    def __init__(self, X, points=100):

        self.points = points

        # calcualte PCA
        self.pca = PCA(X)
        # get kernel for PC1 and PC2
        self.kernel = gaussian_kde(self.pca.T[:,[0,1]].T)

        # calcualte spatial span
        self.X, self.Y = np.mgrid[self.pca.T.min(0)[0]:self.pca.T.max(0)[0]:(self.points*1j),
                                  self.pca.T.min(0)[1]:self.pca.T.max(0)[1]:(self.points*1j)]
        self.positions = np.vstack([self.X.ravel(), self.Y.ravel()])

        self.values = self.kernel(self.positions)

        self.Z = np.reshape(self.values,self.X.shape)

        # one cell area
        self.cell_area = ((self.pca.T[:,:2].max(0) - self.pca.T[:,:2].min(0))/self.points).prod()

    def area(self,fraction=0.9):
        return self.cell_area * (self.values >  self.Z.max()*(1-fraction)).sum()

    def contour(self,
