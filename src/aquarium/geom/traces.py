import numpy as np
from scipy.spatial.distance import cdist, pdist

def diff(trace):
    n= len(trace)
    if n < 2:
        return trace

    trace_diff = []
    for nr,row in enumerate(trace[:-1]):
        trace_diff.append(pdist(np.vstack((row,trace[nr+1])),metric='euclidean'))

    return np.array(trace_diff)




