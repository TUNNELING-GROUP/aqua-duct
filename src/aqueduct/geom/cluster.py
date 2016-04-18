'''
This module provides functions for clusterization.
Clusterization is done by :mod:`scikit-learn` module.
'''

import numpy as np

from sklearn.cluster import DBSCAN, AffinityPropagation, KMeans, MeanShift, estimate_bandwidth
#from sklearn.cluster import MeanShift as _MeanShift

#from sklearn import metrics
#from sklearn.preprocessing import StandardScaler

from aqueduct.utils.helpers import Auto

def MeanShiftBandwidth(X,**kwargs):

    if 'bandwidth' in kwargs:
        if kwargs['bandwidth'] is Auto:
            bandwidth = estimate_bandwidth(np.array(X), quantile=0.5) # this seems to be optimal
            kwargs.update({'bandwidth':bandwidth})
    return kwargs

def perform_clustering(coords,method,**kwargs):
    # special cases
    if method is MeanShift:
        kwargs = MeanShiftBandwidth(coords,**kwargs)
    clust = method(**kwargs).fit(coords)
    return map(int,clust.labels_ + 1) # TODO: add one to it? Is it OK?


