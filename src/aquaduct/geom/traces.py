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

import numpy as np
from scipy.spatial.distance import pdist


################################################################################
# traces utils and helpers

def diff(trace):
    """
    This function calculates the distance between 2 given points.

    :param trace: coordinates in numpy array object
    :return: distance between points
    """
    # trace - 3D coordinates
    # returns distances between coordinates
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
    """
    :param start: coordinates of the first point as a numpy array object
    :param stop: coordinates of the second point as a numpy array object
    :param nr: number of elements between the first and second point
    :return: two-dimentional numpy array; number of dimentions depends on nr parameter
    """
    # returns points between start and stop as linear interpolations
    # if nr == 1 then midpoint is returned
    return np.array([np.linspace(cb, ce, nr + 2)[1:-1] for cb, ce in zip(start, stop)]).T


def midpoints(paths):
    """
    The function returns a tuple of numpy arrays extended with mid point spanning last and first element(column)
     of these arrays.

    :param paths: a tuple of 2-dimentional np.arrays that hold 3D coordinates; each element holds one trace, all elements are supposed to make one path divided in to sections
    :return: paths elements with additional mid points as a generator object
    """
    # paths - a tuple of 2d np.arrays that hold 3D coordinates, each element holds one trace, all elements are supposed to make one path divided in to sections
    # yields paths elements with additional mid points
    # if input paths is follwoing:
    #   11111 33333 55555
    # function yields the same elements plus midpoints:
    #   111112 2333334 455555
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
            next_trace = []
            if nr + 1 < n:
                if paths[nr + 1].size > 0:
                    next_trace = paths[nr + 1][0]
            # calculate midpoints if relevant
            if len(trace) > 0:
                if len(past) > 0:
                    midp = tracepoints(past, trace[0], 1)
                    trace = np.vstack((midp, trace))
                if len(next_trace) > 0:
                    midp = tracepoints(trace[-1], next_trace, 1)
                    trace = np.vstack((trace, midp))
            yield trace
    else:
        yield paths[0]


def length_step_std(trace):
    """
    This function calculates sum, mean and standard deviation from all segments of a trace.

    :param trace: coordinates of points as numpy array
    :return: a tuple with basics statistics of a trace
    """
    # trace - 3D coordinates
    # calculates diff over trace and returns sum, mean and std of diff
    # if trace is empty or have length < 2 nans are returned
    if len(trace) < 2:
        return float('nan'), float('nan'), float('nan')
    d = diff(trace)
    return np.sum(d), np.mean(d), np.std(d)


def derrivative(values):
    # FIXME: This function is porbably completely wrong! AP is going to rewrite it.
    # values - 3D coordinates
    # calculates derrivative of lenght of trace
    # uses diff but yields the same number of values as in input data
    # this is done by interpolation and applying simple correction
    diff = np.diff(values)
    size = len(diff)
    correction = 1. / (size + 1)  # this correct values so after integration they are closer to expected value
    # This was calculated by following experiment:
    # w = []
    # for q in range(10000):
    #     r = np.cumsum(np.random.rand(10000))
    #     w.append(r[-1]/sum(list(traces.derrivative(r))))
    # np.mean(w) # w is of narmal distribution

    for nr in range(size + 1):
        if nr == 0:
            yield diff[0] - correction  # begin
        elif nr == size:
            yield diff[-1] - correction  # end
        else:
            yield (diff[nr - 1] + diff[nr]) / 2. - correction

# todo : aby zaoszczedzic na obliczeniach mozna pomijac takie katy(lub zwracac 0), ktorych zwracane wartosci są bardzo,bardzo małe (rzedu np 10**-4)-> np kat 0.005 rad to 0,29stopnia miary łukowej
# wektory: promień atomu wodoru to 0.529A

################################################################################
# vectors, triangles and angles helpers

def vector_norm(V):
    """

    :param V: a vector in a form of array-like object, tuple or a list
    :return: normalized length of a vector
    """
    # calculate length of physical vector based on it's coordinates
    # input: tuple or a list
    # output: float
    # return np.sqrt(np.dot(V, V.conj()))
    return np.sqrt(np.dot(V, V))
    # return np.linalg.norm(V)


def triangle_angles(A, B, C):
    """
    Parameters are coordinates of points which are tops of triangle. The function calculates angles in a triangle formed by given coordinates.

    :param A: coordinates of the first point
    :param B: coordinates of the second point
    :param C: coordinates of the third point
    :return: list of arguments where angle is given in radians , the output is as follow: [BAC,CAB,ABC]
    """
    # http://stackoverflow.com/questions/5122372/angle-between-points
    # A,B,C are point in the space
    # input: 3 space coords of points (as tuple or list)
    # returns list of arguments where angle is given in radians , the output is as follow: [BAC,CAB,ABC]
    A, B, C = map(np.array, (A, B, C))
    a = C - A
    b = B - A
    c = C - B
    angles = []
    for e1, e2 in ((a, b), (a, c), (b, -c)):
        num = np.dot(e1, e2)
        denom = vector_norm(e1) * vector_norm(e2)
        angles.append(np.arccos(np.clip(num / denom, -1, 1)))  # cliping values in to [-1,1]
        if np.isnan(angles[-1]):
            angles[-1] = 0.
    return angles


def triangle_angles_last(A, B, C):
    """
    Parameters are coordinates of points which are tops of triangle. The function calculates the [ABC] angle.

    :param A: coordinates of the first point [A top]
    :param B: coordinates of the second point [B top]
    :param C: coordinates of the third point [C top]
    :return: list with one value of ABC angle in radians
    """
    # http://stackoverflow.com/questions/5122372/angle-between-points
    # A,B,C are point in the space
    # input: 3 space coords of points (as tuple or list)
    # returns list with one value of ABC angle in radians
    A, B, C = map(np.array, (A, B, C))
    a = C - A
    b = B - A
    c = C - B
    angles = []
    for e1, e2 in ((b, -c),):
        num = np.dot(e1, e2)
        denom = vector_norm(e1) * vector_norm(e2)
        angles.append(np.arccos(np.clip(num / denom, -1, 1)))  # cliping values in to [-1,1]
    return angles


def triangle_height(A, B, C):
    """
    Parameters are coordinates of points which are tops of triangle. The function calculates the ABC triangle height.

    :param A: coordinates of the first point [A top]
    :param B: coordinates of the second point [B top]
    :param C: coordinates of the third point [C top]
    :return: one value of ABC triangle height
    """
    # a is head
    # input: 3 space coords of points (as tuple or list)
    # output float, value of triangle height
    angles = triangle_angles_last(A, B, C)  # ta funkcja zwraca tylko 1 element
    A, B, C = map(np.array, (A, B, C))
    c = vector_norm(B - A)
    h = np.sin(angles[-1]) * c
    if np.isnan(h):
        h = 0.
    return h


def vectors_angle(A, B):
    """
    This function calculates the angle between two given vectors (starting from the [0,0,0] to the given coordinates.

    :param A: coordinates of the first point which is the end of the vector
    :param B: coordinates of the second point which is the end of the vector
    :return: the angle between vectors in question (in radians)
    """
    angle = np.arccos(np.dot(A, B) / (vector_norm(A) * vector_norm(B)))
    if np.isnan(angle):
        return 0.0
    return angle


def vectors_angle_alt(A, B):
    '''
    This function calculates the angle between two given vectors (starting from the [0,0,0] to the given coordinates

    * alternative method.

    :param A: coordinates of the first point which is the end of the vector
    :param B: coordinates of the second point which is the end of the vector
    :return: the angle between vectors in question (in radians)
    '''
    return np.arccos(np.clip(np.dot(A / vector_norm(A), B / vector_norm(B)), -1.0, 1.0))


def vectors_angle_alt_anorm(A, B, A_norm):
    """
    This function calculates the angle between two given vectors (starting from the [0,0,0] to the given coordinates
     - alternative method with additional A_norm holding norm of A.

    :param A: coordinates of the first point which is the end of the vector
    :param B: coordinates of the second point which is the end of the vector
    :param A_norm: additional parameter holding normalized of vector A
    :return: the angle between vectors in question (in radians)
    """
    return np.arccos(np.clip(np.dot(A / A_norm, B / vector_norm(B)), -1.0, 1.0))


def vectors_angle_anorm(A, B, A_norm):
    """
    This function calculates the angle between two given vectors (starting from the [0,0,0] to the given coordinates
     using additional A_norm holding norm of A.

    :param A: coordinates of the first point which is the end of the vector
    :param B: coordinates of the second point which is the end of the vector
    :param A_norm: additional parameter holding normalized of vector A
    :return: the angle between vectors in question (in radians)
    """
    norm2 = A_norm * vector_norm(B)
    if norm2 == 0.:
        return 0.
    angle = np.clip(np.dot(A, B) / norm2, -1., 1.)
    if np.isnan(angle):
        return 0.
    return np.arccos(angle)


################################################################################
# Linearization classes

# poprawiony algorytm
class LinearizeOneWay(object):
    def here(self, coords):
        # coords - 3D coordintates of a trace
        # yields indices of coords which is a staring point of linear fragments of the trace and the next point after
        # enf of linear segment; done in one way
        """
        This function simplifies the trace by removing the redundant, linear points
        :param coords: 3D coordinates of a trace as an array-like object
        :return: indices of coordinates which are a staring and ending points of linear fragments and other non-linear
        points of the trace
        """
        size = len(coords)
        yield 0
        if size <= 3:
            yield 2
        else:
            for sp in range(size - 2):
                ep = sp + 3
                if self.is_linear(coords[sp:ep]):
                    continue
                yield sp + 1
                continue
            yield size - 1


class LinearizeHobbit(LinearizeOneWay):
    def and_back_again(self, coords):
        # coords - 3D coordintates of a trace
        # yields indices of coords that spans linear fragments of the trace; done in opposite way than in one way
        size = len(coords)
        return (size - e - 1 for e in self.here(coords[::-1]))

    def __call__(self, coords):
        # coords - 3D coordintates of a trace
        # wrapper that uses here and and_back_again methods to get merged uniq and sorted indices of coords that spans linear fragments of the trace
        # returns these points from coords that are linear simplification of coords
        # __call__ is required by child classes
        here = self.here(coords)
        and_back_again = self.and_back_again(coords)
        linearize = sorted(list(set(list(here) + list(and_back_again))))
        return coords[linearize]


class LinearizeRecursive(object):
    """
    Base class for linearization methods classes.

    It implements recursive algorithm.
    """

    def here(self, coords, depth=0):
        """
        Core of recursive linearization argorithm.

        It checks if the first, the last and the middle point are linear according to the criterion. The middle point is a selected point that is in the middle of length of the paths made by input coordinates.

        If these points are linear their indices are returned. Otherwise, coordinates are split into two parts. First part spans points from the first point to the middle point (inclusive) and the second part spans points from the middle (inclusive) to the last point. Next, these two parts are submitted recursively to :meth:`here`.

        Results of these recursive calls are joined, redundant indices are removed and sorted result is returned.

        :param numpy.ndarray coords: Input coordinates.
        :param int depth: Depth of recurence.
        :return: Indices of ``coords`` points that can be used instead of all points in visulatization.
        :rtype: list of int
        """
        # klasa nie ma zdefiniowanej metody is_linear
        depth += 1

        lengths = np.hstack(([0], np.cumsum(diff(coords))))
        size = len(lengths)
        if size <= 3:
            return range(size)
        sp = 0
        ep = size - 1
        mp = np.argwhere(lengths > max(lengths) / 2)
        if len(mp): # in case lenghs are zeros
            mp = int(mp[0])
        else:
            mp = (sp + ep) / 2

        if mp == sp:
            mp += 1
        if mp == ep:
            mp -= 1

        if self.is_linear(coords[[sp, mp, ep]], depth=depth):
            return [sp, mp, ep]
        return sorted(
            list(set(self.here(coords[:mp + 1], depth=depth) + [e + mp for e in self.here(coords[mp:], depth=depth)])))

    def __call__(self, coords):
        # returns these points from coords that are linear simplification of coords
        # __call__ is required by child classes
        here = self.here(coords)
        return coords[here]


################################################################################
# Linearize criteria

class TriangleLinearize(object):
    def __init__(self, threshold=0.01):
        # threshold - maximal allowed sum of heights of triangles made of beginning, end and all middle points
        ## bardzo ostre kryterium!!
        self.threshold = threshold

    def is_linear(self, coords, **kwargs):
        # coords - 3D coordintates of a trace
        # returns True if coords make a straight line
        # criterion of linearity:
        # if sum of heights of triangles made of beginning, end and all middle points does not exceed threshold coords are linear
        list_of_h = list()
        for head in coords[1:-1]:
            list_of_h.append(triangle_height(head, coords[0], coords[-1]))
            # print list_of_h, sum(list_of_h)
            if sum(list_of_h) > self.threshold:
                return False
        return True


class VectorLinearize(object):
    """
    Base class for linearization methods classes.

    It implements vector linearization criterion.
    """

    def __init__(self, treshold=0.05236):
        # TODO: Add docs
        self.treshold = treshold

    def is_linear_core(self, coords, depth=None):
        """
        Method checks if input coordinates are linear according to the threshold and depth.

        It begins with calculation of the threshold. If `depth` is None it is set to 1. Current threshold is calculated with following simple equation:

        .. math::

            threshold_{current} = threshold_{initial} * (2 - 0.9^{depth})

        Next, in a loop over all points but the first and the last the angle is calculated between two vectors. The first one made by the point and the first point, and the second vector made by the last and the first point. If any of the calculated angles is bigger the the treshold methods returns False; otherwise method returns True.

        :param numpy.ndarray coords: Coordinates for which linearization criterion is checked.
        :param int depth: Depth of recurence.
        :return: True if input coordinates are linear and False otherwise.
        :rtype: bool
        """
        if depth is None:
            depth = 1

        treshold = self.treshold + self.treshold * (1 - 0.9 ** depth)  # FIXME: magic constant!

        V = coords[-1] - coords[0]
        V_norm = vector_norm(V)
        for cp in coords[:-1]:
            V_sum = cp - coords[0]
            if vectors_angle_anorm(V, V_sum, V_norm) > treshold:
                return False
        return True

    def is_linear(self, coords, depth=None, **kwargs):

        """
        For more detail see :meth:`is_linear_core` which is used as the criterion of linearity in this method.

        :param numpy.ndarray coords: Coordinates for which linearization criterion is checked.
        :param int depth: Depth of recurence.
        :return: True if input coordinates are linear and False otherwise. Criterion is checked for coordinates in normal and reverse order.
        :rtype: bool
        """
        if not self.is_linear_core(coords, depth=depth):
            return False
        elif not self.is_linear_core(coords[::-1], depth=depth):
            return False
        return True


################################################################################
# Concrete classes for linearization

class LinearizeRecursiveVector(LinearizeRecursive, VectorLinearize):
    """
    ..  _simply_smooths_details:

    Class provides recursive linearization of coordinates with :class:`LinearizeRecursive` algorithm and the criterion of linearity implemented by :class:`VectorLinearize`. This is default method.
    """
    pass


class LinearizeRecursiveTriangle(LinearizeRecursive, TriangleLinearize):
    """
    Class provides recursive linearization of coordinates with :class:`LinearizeRecursive` algorithm and the criterion of linearity implemented by :class:`TriangleLinearize`.
    """
    pass


class LinearizeHobbitVector(LinearizeHobbit, VectorLinearize):
    """
    Class provides recursive linearization of coordinates with :class:`LinearizeHobbit` algorithm and the criterion of linearity implemented by :class:`VectorLinearize`.
    """
    pass


class LinearizeHobbitTriangle(LinearizeHobbit, TriangleLinearize):
    """
    Class provides recursive linearization of coordinates with :class:`LinearizeHobbit` algorithm and the criterion of linearity implemented by :class:`TriangleLinearize`.
    """
    pass


class LinearizeOneWayVector(LinearizeOneWay, VectorLinearize):
    """
    Class provides recursive linearization of coordinates with :class:`LinearizeOneWay` algorithm and the criterion of linearity implemented by :class:`VectorLinearize`.
    """
    pass


class LinearizeOneWayTriangle(LinearizeOneWay, TriangleLinearize):
    """
    Class provides recursive linearization of coordinates with :class:`LinearizeOneWay` algorithm and the criterion of linearity implemented by :class:`TriangleLinearize`.
    """
    pass

################################################################################
