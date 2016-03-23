'''
This module provides functions for clusterization.
Clusterization is done by :mod:`scikit-learn` module.
'''

import numpy as np

from sklearn.cluster import DBSCAN
#from sklearn import metrics
#from sklearn.preprocessing import StandardScaler



# lets use DBSCAN as default method

def get_default_culstering_method(**kwargs):
    '''
    this is very simple
    '''
    return DBSCAN(**kwargs)

def perform_clustering(coords,**kwargs):

    clust = get_default_culstering_method(**kwargs).fit(coords)

    return clust.labels_


