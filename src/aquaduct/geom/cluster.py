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

"""
This module provides functions for clusterization.
Clusterization is done by :mod:`scikit-learn` module.
"""

import numpy as np
from sklearn.cluster import Birch, DBSCAN, AffinityPropagation, KMeans, MeanShift, estimate_bandwidth
from scipy.spatial.distance import cdist

# problems with clustering methods and size of set
# DBSCAN:              n > 0
# AffinityPropagation: n > 0
# KMeans:              n > clusters
# MeanShift:           n > 6


AVAILABLE_METHODS = ['dbscan','kmeans', 'affprop', 'meanshift', 'birch','barber']

def get_required_params(method):
    if method == 'kmeans':
        return ['n_clusters']


from aquaduct.utils.helpers import Auto
from aquaduct.utils import clui
from aquaduct.traj.barber import WhereToCut

class BarberClusterResult(object):
    '''
    Helper class for results of barber clusterization.
    '''
    def __init__(self,labels_):
        self.labels_ = np.array(labels_)

class BarberCluster(object):
    '''
    Wrapper class that implements *barber* clusterization.
    '''
    def fit(self, coords, spheres=None):
        '''
        :param Iterable coords: Input coordinates of points to be clustered.
        :param Iterable spheres: Input spheres for each point.
        '''
        wtc = WhereToCut()
        wtc.spheres = spheres
        clouds = wtc.cloud_groups(progress=True)
        # clouds to labels!
        labels = np.zeros(len(spheres))
        spheres_id = [s.nr for s in spheres]
        for cloud_id,cloud in clouds.iteritems():
            labels[[spheres_id.index(c) for c in cloud]] = cloud_id
        return BarberClusterResult(labels)


def MeanShiftBandwidth(X, **kwargs):
    '''
    Helper function for automatic calculation of a bandwidth for MeanShift method.

    :param Iterable X: Coordinates of points to be clustered.
    '''
    if 'bandwidth' in kwargs:
        if kwargs['bandwidth'] is Auto:
            bandwidth = estimate_bandwidth(np.array(X), quantile=0.5)  # TODO: change it to the default value of 0.3 or use it as option?
            kwargs.update({'bandwidth': bandwidth})
            clui.message("Meanshift automatic bandwidth calculation: bandwidth = %f" % float(
                bandwidth))  # TODO: make it properly
    return kwargs


class PerformClustering(object):
    '''
    Helper class for clusterization.
    '''
    # aqeuduct clustering helper class

    def __init__(self, method, **kwargs):
        '''
        :param object method: Class that implements cclusterization via *fit* method.
        '''

        self.method = method
        self.method_kwargs = kwargs
        self.method_results = None
        self.clusters = None

    def __str__(self):
        out = str(self.method.__name__)
        return out

    def __call__(self, coords, spheres=None):
        # compatibility
        return self.fit(coords, spheres=spheres)

    def _get_noclusters(self, n):
        return [0] * n

    def fit(self, coords, spheres=None):
        '''
        :param Iterable coords: Input coordinates of points to be clustered.
        :param Iterable spheres: Input spheres for each point. Optional, important only if :attr:`method` is :class:`BarberCluster`.
        :return: Clusters numbers.
        :rtype: list of int
        '''
        # spheres are used for Barber only
        if len(coords) < 2:
            self.clusters = self._get_noclusters(len(coords))
            return self.clusters
        # special cases
        if self.method is BarberCluster:
            method = self.method()
            self.method_results = method.fit(coords, spheres)
            self.clusters = map(int, self.method_results.labels_ + 1)
        else:
            if self.method is MeanShift:
                if len(coords) < 6:
                    self.clusters = self._get_noclusters(len(coords))
                    clui.message("Number of objects %d < 6 is too small for MeanShift method. Skipping." % (len(coords)))
                    return self.clusters
                method = self.method(**MeanShiftBandwidth(coords, **self.method_kwargs))
            else:
                if self.method is KMeans:
                    if 'n_clusters' in self.method_kwargs:
                        if len(coords) < self.method_kwargs['n_clusters']:
                            self.clusters = self._get_noclusters(len(coords))
                            clui.message(
                                "Number of objects %d < %d is too small for KMeans method. Skipping." % (
                                    len(coords), self.method_kwargs['n_clusters']))
                            return self.clusters
                method = self.method(**self.method_kwargs)
            self.method_results = method.fit(coords)
            self.clusters = map(int, self.method_results.labels_ + 1)
        return self.clusters

    def centers(self):
        '''
        :return: Centers of clusters.
        '''
        assert self.clusters is not None, "Perform clusterization first."
        assert self.method_results is not None, "Perform clusterization first."

        if hasattr(self.method_results, 'cluster_centers_'):
            return self.method_results.cluster_centers_
        raise NotImplementedError('Cluster centers is not implemented for %r method yet.' % self.method)
