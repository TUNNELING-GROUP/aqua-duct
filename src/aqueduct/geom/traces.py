import numpy as np
from scipy.spatial.distance import cdist, pdist


def diff(trace):
    assert isinstance(trace, np.ndarray), "Trace should be of np.ndarray type, %r submited instead." % type(trace)
    assert len(trace.shape) == 2, "Traces should be 2d, %dd submited instead (trace no. %d)" % len(trace.shape)

    n = len(trace)
    if n < 2:
        return None

    trace_diff = []
    for nr, row in enumerate(trace[:-1]):
        trace_diff.append(float(pdist(np.vstack((row, trace[nr + 1])), metric='euclidean')))

    return np.array(trace_diff)


def tracepoints(start, stop, nr):
    # nr == 1 then midpoint is returned
    return np.array([np.linspace(cb, ce, nr + 2)[1:-1] for cb, ce in zip(start, stop)]).T


def midpoints(paths):
    # paths is a tuple of 2d np.arrays,
    assert isinstance(paths, tuple), "Paths should be of tuple type, %r submitted instead." % type(paths)
    for nr, trace in enumerate(paths):
        assert isinstance(trace,
                          np.ndarray), "Traces should be of numpy.ndarray type, %r submited instead (trace no. %d)" % (
            type(trace), nr)
        # assert len(trace.shape) == 2, "Traces should be 2d, %dd submited instead (trace no. %d)" %(len(trace.shape),nr)

    n = len(paths)
    last_trace = None
    if n > 1:
        for nr, trace in enumerate(paths):
            # find past and next trace
            past = []
            if nr - 1 >= 0:
                if paths[nr - 1].size > 0:
                    past = paths[nr - 1][-1]
            next = []
            if nr + 1 < n:
                if paths[nr + 1].size > 0:
                    next = paths[nr + 1][0]
            # calculate midpoints if relevant
            if len(trace) > 0:
                if len(past) > 0:
                    midp = tracepoints(past, trace[0], 1)
                    trace = np.vstack((midp, trace))
                if len(next) > 0:
                    midp = tracepoints(trace[-1], next, 1)
                    trace = np.vstack((trace, midp))
            yield trace
    else:
        yield paths[0]


def length_step_std(trace):
    d = diff(trace)
    return np.sum(d),np.mean(d),np.std(d)
