import numpy as np
from scipy.spatial.distance import cdist, pdist

def diff(trace):
    n = len(trace)
    if n < 2:
        return None

    trace_diff = []
    for nr,row in enumerate(trace[:-1]):
        trace_diff.append(float(pdist(np.vstack((row,trace[nr+1])),metric='euclidean')))

    return np.array(trace_diff)

def midpoints(path):
    n = len(path)
    for nr,trace in enumerate(path):
        # mid points!
        if len(trace) > 0:
            if nr == 0:
                if len(path[nr+1]) > 0:
                    midp = np.mean(np.vstack((trace[-1], path[nr + 1][0])), 0)
                    trace = np.vstack((trace,midp))
            elif nr == n-1:
                if len(path[nr-1]) > 0:
                    midp = np.mean(np.vstack((trace[0], path[nr - 1][-1])), 0)
                    trace = np.vstack((midp,trace))
            else:
                if len(path[nr-1]) > 0:
                    midp = np.mean(np.vstack((trace[0], path[nr - 1][-1])), 0)
                    trace = np.vstack((midp,trace))
                if len(path[nr+1]) > 0:
                    midp = np.mean(np.vstack((trace[-1], path[nr + 1][0])), 0)
                    trace = np.vstack((trace,midp))
        yield trace


