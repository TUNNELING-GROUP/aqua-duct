import numpy as np
from scipy.spatial.distance import cdist, pdist
import copy

from aqueduct.utils.helpers import arrayify,lind


def vector_norm(V):
    # return np.sqrt(np.dot(V, V.conj()))
    return np.sqrt(np.dot(V, V))
    # return np.linalg.norm(V)


def triangle_angles(A, B, C):
    # http://stackoverflow.com/questions/5122372/angle-between-points
    # ABC are point in the space
    A, B, C = map(np.array, (A, B, C))
    a = C - A
    b = B - A
    c = C - B
    angles = []
    for e1, e2 in ((a, b), (a, c), (b, -c)):
        num = np.dot(e1, e2)
        denom = vector_norm(e1) * vector_norm(e2)
        angles.append(np.arccos(num / denom))
        if np.isnan(angles[-1]):
            angles[-1] = 0.
    return angles


def triangle_angles_last(A, B, C):
    # http://stackoverflow.com/questions/5122372/angle-between-points
    # ABC are point in the space
    A, B, C = map(np.array, (A, B, C))
    a = C - A
    b = B - A
    c = C - B
    angles = []
    for e1, e2 in ((b, -c),):
        num = np.dot(e1, e2)
        denom = vector_norm(e1) * vector_norm(e2)
        angles.append(np.arccos(num / denom))
    return angles


def triangle_height(A, B, C):
    # a is head
    angles = triangle_angles_last(A, B, C)
    A, B, C = map(np.array, (A, B, C))
    c = vector_norm(B - A)
    h = np.sin(angles[-1]) * c
    if np.isnan(h):
        h = 0.
    return h


def vectors_angle(A, B):
    angle = np.arccos(np.dot(A, B) / (vector_norm(A) * vector_norm(B)))
    if np.isnan(angle):
        return 0.0
    return angle


def vectors_angle_alt(A, B):
    return np.arccos(np.clip(np.dot(A / vector_norm(A), B / vector_norm(B)), -1.0, 1.0))


def vectors_angle_alt_anorm(A, B, A_norm):
    return np.arccos(np.clip(np.dot(A / A_norm, B / vector_norm(B)), -1.0, 1.0))


def vectors_angle_anorm(A, B, A_norm):
    angle = np.arccos(np.dot(A, B) / (A_norm * vector_norm(B)))
    if np.isnan(angle):
        return 0.0
    return angle



class LinearizeOneWay(object):

    def here(self, coords):
        size = len(coords)
        yield 0
        ep = 0
        for sp in range(size):
            if sp < ep:
                continue
            for ep in range(sp + 2, size):
                if self.is_linear(coords[sp:ep]):
                    continue
                yield ep
                break


class LinearizeHobbit(LinearizeOneWay):


    def and_back_again(self, coords):
        size = len(coords)
        return (size - e - 1 for e in self.here(coords[::-1]))

    def __call__(self, coords):
        here = self.here(coords)
        and_back_again = self.and_back_again(coords)
        linearize = sorted(list(set(list(here) + list(and_back_again))))
        return coords[linearize]

class LinearizeRecursive(object):

    def here(self,coords,depth=0):
        depth += 1

        lengths = np.hstack(([0],np.cumsum(diff(coords))))
        size = len(lengths)
        if size <= 3:
            return range(size)
        sp = 0
        ep = size - 1
        mp = int(np.argwhere(lengths>max(lengths)/2)[0])

        if mp == sp:
            mp += 1
        if mp == ep:
            mp -= 1

        if self.is_linear(coords[[sp,mp,ep]]):
            return [sp,mp,ep]
        return sorted(list(set(self.here(coords[sp:mp+1],depth=depth) + [e+mp for e in self.here(coords[mp:ep],depth=depth)])))


    def __call__(self, coords):
        here = self.here(coords)
        return coords[here]


class TrianlgeLinearize(object):
    def __init__(self, treshold):

        self.treshold = treshold

    def is_linear(self, coords):
        list_of_h = list()
        for head in coords[1:-1]:
            list_of_h.append(triangle_height(head, coords[0], coords[-1]))
            # print list_of_h, sum(list_of_h)
            if sum(list_of_h) > self.treshold:
                return False
        return True


class VectorLinearize(object):
    def __init__(self, treshold):

        self.treshold = treshold

    def is_linear_core(self,coords):
        V = coords[-1] - coords[0]
        V_norm = vector_norm(V)
        for cp in coords[:-1]:
            V_sum = cp - coords[0]
            if vectors_angle_anorm(V, V_sum, V_norm) > self.treshold:
                return False
        return True

    def is_linear(self, coords):
        if not self.is_linear_core(coords):
            return False
        elif not self.is_linear_core(coords[::-1]):
            return False
        return True


class LinearizeRecursiveVector(LinearizeRecursive,VectorLinearize):
    pass




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
    return np.sum(d), np.mean(d), np.std(d)
