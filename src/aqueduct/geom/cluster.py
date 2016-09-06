# -*- coding: utf-8 -*-
'''
This module provides functions for clusterization.
Clusterization is done by :mod:`scikit-learn` module.
'''

import numpy as np
from sklearn.cluster import Birch, DBSCAN, AffinityPropagation, KMeans, MeanShift, estimate_bandwidth

# problems with clustering methods and size of set
# DBSCAN:              n > 0
# AffinityPropagation: n > 0
# KMeans:              n > clusters
# MeanShift:           n > 6


from aqueduct.utils.helpers import Auto
from aqueduct.utils import log

def MeanShiftBandwidth(X, **kwargs):
    if 'bandwidth' in kwargs:
        if kwargs['bandwidth'] is Auto:
            bandwidth = estimate_bandwidth(np.array(X), quantile=0.5)  # this seems to be optimal
            kwargs.update({'bandwidth': bandwidth})
            log.message("\nMeanshift automatic bandwidth calculation: bandwidth = %f" % float(bandwidth)) # TODO: make it properly
    return kwargs



class PerformClustering(object):
    # aqeuduct clustering helper class

    def __init__(self, method, **kwargs):

        self.method = method
        self.method_kwargs = kwargs
        self.method_results = None
        self.clusters = None

    def __call__(self, coords):
        # compatibility
        return self.fit(coords)

    def _get_noclusters(self,n):
        return [0]*n

    def fit(self, coords):
        if len(coords) < 2:
            self.clusters = self._get_noclusters(len(coords))
            return self.clusters
        # special cases
        if self.method is MeanShift:
            if len(coords) < 6:
                self.clusters = self._get_noclusters(len(coords))
                return self.clusters
            method = self.method(**MeanShiftBandwidth(coords, **self.method_kwargs))
        else:
            if self.method is KMeans:
                if 'n_clusters' in self.method_kwargs:
                    if len(coords) < self.method_kwargs['n_clusters']:
                        self.clusters = self._get_noclusters(len(coords))
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
