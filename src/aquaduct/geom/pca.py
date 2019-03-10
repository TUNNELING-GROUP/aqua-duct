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
import scipy.linalg


class NullPrepocess(object):
    def __init__(self):
        pass

    def build(self, X):
        pass

    def __call__(self, X):
        return X

    def undo(self, X):
        return X


class Center(object):
    def __init__(self):
        pass

    def build(self, X):
        self.mean = np.mean(X, 0)

    def __call__(self, X):
        return X - self.mean

    def undo(self, X):
        return X + self.mean


class Normalize(object):
    def __init__(self):
        pass

    def build(self, X):
        self.std = np.std(X, 0)

    def __call__(self, X):
        return X / self.std

    def undo(self, X):
        return X * self.std


class Standartize(object):
    def __init__(self):
        self.center = Center()
        self.normalize = Normalize()

        test_arr = np.random.random((10, 3))
        test_cntr = np.random.random((1, 3))
        P = Polarize(center=test_cntr)
        P.build(test_arr)
        test_arr_p = P(test_arr)
        Pp = P.undo(test_arr_p) - test_arr
        [self.assertAlmostEqual(x, 0, 7) for x in Pp.ravel()]

    def test_polarize_in(self):
        P = Polarize(center=np.random.random((1, 3)))
        self.assertRaises(TypeError, P, 'cupkaces')
        self.assertRaises(TypeError, Polarize, 'cupcakes')
        self.assertRaises(TypeError, P, np.random.random((1, 2)))


if __name__ == '__main__':
    def build(self, X):
        self.center.build(X)
        self.normalize.build(X)


    def __call__(self, X):
        return self.normalize(self.center(X))


    def undo(self, X):
        return self.center.undo(self.normalize.undo(X))


class Polarize(object):
    # TODO: equal variance of t and f
    def __init__(self, center=np.array([0, 0, 0]), rvar=0.1, equaltf=True):
        '''
        Prepocessing filter for 3D cartesian coordinates transformation to spherical coordinates.

        .. note::

            Component f is in range 0 - 2\pi.

        :param center: Center of the hypothetical sphere.
        :param rvar: Desired amount of variance of *r* component measured as fraction of mean *t* and *f* variance.
        :param equaltf: If set ``True``, *t* range is scaled to *f*.
        '''

        if not isinstance(center, np.ndarray):
            raise TypeError('Constructor called with center param of invalid type')
        elif np.shape(center) != (3,):
            raise TypeError('Constructor called with center param of invalid shape')
        else:
            self.center = center

        self.rvar = rvar
        self.equaltf = equaltf
        self.tmean = 0
        self.fmean = 0
        self.rvar_factor = 1

    @property
    def _mt(self):
        return 2. if self.equaltf else 1.

    def _Xrtf(self, X):
        X = X - self.center
        # ISO, f + pi
        r = (X ** 2).sum(1) ** 0.5
        t = np.arccos(X[:, 2] / r)  # 0,pi
        f = np.arctan2(X[:, 1], X[:, 0]) + np.pi  # -pi,pi # f+pi 0,2pi
        return X, r, t, f

    def _circle_tf(self, t, f):
        return t % np.pi, f % (2 * np.pi)

    def build(self, X):
        X, r, t, f = self._Xrtf(X)
        # calculate mean values for polar components (circle it)
        self.tmean, self.fmean = self._circle_tf(np.pi / 2 - np.mean(t), np.pi - np.mean(f))
        tf_var = (np.var(t) * 2 + np.var(f)) / 4.
        if np.std(r) > 0:
            self.rvar_factor = ((tf_var * self.rvar) ** 0.5) * (1. / np.std(r))
        else:
            self.rvar_factor = 1.

    def __call__(self, X):
        X, r, t, f = self._Xrtf(X)
        # center t and f with appropriate mean values
        t += self.tmean
        f += self.fmean
        # take care of polarity
        t, f = self._circle_tf(t, f)
        return np.array([r * self.rvar_factor, t * self._mt, f]).T

    def undo(self, X):
        # remove centering
        t = X[:, 1] / self._mt - self.tmean
        f = X[:, 2] - self.fmean
        # take care of polarity
        t, f = self._circle_tf(t, f)
        # correct f to iso
        f -= np.pi
        # get xyz
        x = X[:, 0] / self.rvar_factor * np.sin(t) * np.cos(f)
        y = X[:, 0] / self.rvar_factor * np.sin(t) * np.sin(f)
        z = X[:, 0] / self.rvar_factor * np.cos(t)
        return np.array([x, y, z]).T + self.center


class PCA(object):
    def __init__(self, preprocess=NullPrepocess()):
        self.preprocess = preprocess

    def build(self, X):
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
            return np.dot(self.preprocess(X), self.P[:, pc])
        return np.dot(self.preprocess(X), self.P)

    def undo(self, T, pc=None):
        if pc:
            return self.preprocess.undo(np.dot(T, self.Pt[pc, :]))
        return self.preprocess.undo(np.dot(T, self.Pt))
