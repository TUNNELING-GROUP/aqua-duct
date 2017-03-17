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

class BarberClusterResult(object):

    def __init__(self,labels_):
        self.labels_ = np.array(labels_)

class BarberCluster(object):

    def fit(self,coords,radii=None):
        friends = {}
        for nr,(coord,radius) in enumerate(zip(coords,radii)):
            distances = cdist([coord],coords,metric='euclidean').flatten() - np.array(radii) - radius
            # less then zero are intersecting
            if (distances<=0).any():
                friends.update({nr:np.argwhere(distances<=0).flatten().tolist()})
        # loop over friends' groups
        clustered = []
        clusters = []
        for leader,group in friends.iteritems():
            if leader in clustered: continue
            cluster = [] + group
            last_size = -1
            while len(cluster) != last_size:
                last_size = len(cluster)
                for leader2,group2 in friends.iteritems():
                    if leader2 in clustered: continue
                    if set(cluster).intersection(set(group2)):
                        cluster += group2
                        cluster = list(set(cluster))
            clustered += cluster
            clusters.append(cluster)
        clusters.sort(key=lambda x: len(x),reverse=True)
        # make labels
        labels = [None]*len(coords)
        for nr,cluster in enumerate(clusters):
            for c in cluster:
                labels[c] = nr
        return BarberClusterResult(labels)


def MeanShiftBandwidth(X, **kwargs):
    if 'bandwidth' in kwargs:
        if kwargs['bandwidth'] is Auto:
            bandwidth = estimate_bandwidth(np.array(X), quantile=0.5)  # TODO: change it to the default value of 0.3 or use it as option?
            kwargs.update({'bandwidth': bandwidth})
            clui.message("Meanshift automatic bandwidth calculation: bandwidth = %f" % float(
                bandwidth))  # TODO: make it properly
    return kwargs


class PerformClustering(object):
    # aqeuduct clustering helper class

    def __init__(self, method, **kwargs):

        self.method = method
        self.method_kwargs = kwargs
        self.method_results = None
        self.clusters = None

    def __str__(self):
        out = str(self.method.__name__)
        return out

    def __call__(self, coords, radii=None):
        # compatibility
        return self.fit(coords, radii=radii)

    def _get_noclusters(self, n):
        return [0] * n

    def fit(self, coords, radii=None):
        # radii are used for Barber only
        if len(coords) < 2:
            self.clusters = self._get_noclusters(len(coords))
            return self.clusters
        # special cases
        if self.method is BarberCluster:
            method = self.method()
            self.method_results = method.fit(coords,radii)
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
        assert self.clusters is not None, "Perform clusterization first."
        assert self.method_results is not None, "Perform clusterization first."

        if hasattr(self.method_results, 'cluster_centers_'):
            return self.method_results.cluster_centers_
        raise NotImplementedError('Cluster centers is not implemented for %r method yet.' % self.method)

