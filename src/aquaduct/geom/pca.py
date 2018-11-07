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
import scipy.linalg

class NullPrepocess(object):
    def __init__(self):
        pass

    def build(self,X):
        pass

    def __call__(self, X):
        return X

    def undo(self, X):
        return X


class Center(object):
    def __init__(self):
        pass

    def build(self,X):
        self.mean = np.mean(X, 0)

    def __call__(self, X):
        return X - self.mean

    def undo(self, X):
        return X + self.mean


class Normalize(object):
    def __init__(self):
        pass

    def build(self,X):
        self.std = np.std(X, 0)

    def __call__(self, X):
        return X / self.std

    def undo(self, X):
        return X * self.std


class Standartize(object):
    def __init__(self):
        self.center = Center()
        self.normalize = Normalize()

    def build(self,X):
        self.center.build(X)
        self.normalize.build(X)

    def __call__(self, X):
        return self.normalize(self.center(X))

    def undo(self, X):
        return self.center.undo(self.normalize.undo(X))

class Polarize(object):
    def __init__(self,center=np.array([0,0,0]),rvar=0.1):
        self.center = center
        self.rvar = rvar

    def _Xrtf(self,X):
        X = X - self.center
        r = (X**2).sum(1)**0.5
        return X,r,np.arccos(X[:,2]/r),np.arctan2(X[:,1],X[:,0]) + np.pi

    def _circle_tf(self,t,f):
        return t % np.pi, f % (2*np.pi)

    def build(self,X):
        X,r,t,f = self._Xrtf(X)
        # calculate mean values for polar components
        self.tmean = (np.pi/2 - np.mean(t)) % np.pi
        self.fmean = (np.pi - np.mean(f)) % (2*np.pi)
        tf_var = np.var(t)+np.var(f)
        self.rvar_factor = ((tf_var*self.rvar)**0.5)*(1./np.std(r))

    def __call__(self, X):
        X,r,t,f = self._Xrtf(X)
        # center t and f with appropriate mean values
        t += self.tmean
        f += self.fmean
        # take care of polarity
        t,f = self._circle_tf(t,f)
        return np.array([r*self.rvar_factor,t,f]).T

    def undo(self, X):
        # remove centering
        t = X[:,1] - self.tmean
        f = X[:,2] - self.fmean
        # take care of polarity
        t,f = self._circle_tf(t,f)
        # correct f to iso
        f -= np.pi
        # get xyz
        x = X[:,0]/self.rvar_factor*np.sin(t)*np.cos(f)
        y = X[:,0]/self.rvar_factor*np.sin(t)*np.sin(f)
        z = X[:,0]/self.rvar_factor*np.cos(t)
        return np.array([x,y,z]).T + self.center


class PCA(object):
    def __init__(self, preprocess=NullPrepocess()):
        self.preprocess = preprocess

    def build(self,X):
        self.preprocess.build(X)
        # SVD
        self.U, self.d, self.Pt = scipy.linalg.svd(self.preprocess(X), full_matrices=False)
        assert np.all(self.d[:-1] >= self.d[1:]), "SVD error."  # sorted
        # PCA
        self.T = self.U * self.d
        self.eigen = self.d ** 2
        self.sumvariance = np.cumsum(self.eigen)
        self.sumvariance /= self.sumvariance[-1]

    @property
    def P(self):
        return self.Pt.T

    def __call__(self, X, pc=None):
        if pc:
            return np.dot(self.preprocess(X), self.P[:,pc])
        return np.dot(self.preprocess(X), self.P)

    def undo(self, T, pc=None):
        if pc:
            return self.preprocess.undo(np.dot(T, self.Pt[pc,:]))
        return self.preprocess.undo(np.dot(T, self.Pt))
