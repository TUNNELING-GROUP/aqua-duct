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

'''
Smooth module defines methods for smoothing of trajectories.

Available methods:

.. autosummary::
    :nosignatures:

    SavgolSmooth
    WindowSmooth
    DistanceWindowSmooth
    ActiveWindowSmooth
    MaxStepSmooth
    WindowOverMaxStepSmooth
    DistanceWindowOverMaxStepSmooth
    ActiveWindowOverMaxStepSmooth
'''

import numpy as np
from scipy.signal import savgol_filter
from functools import partial

from aquaduct.utils.maths import make_default_array
from aquaduct.geom import traces
from aquaduct.utils.helpers import arrayify


class Smooth(object):
    '''
    Base class for all smoothing methods.
    '''
    def __init__(self, recursive=None, **kwargs):
        '''
        :param int recursive: Number of recursions of the method, everything evaluated to ``False`` is equivalent to 1.
        '''
        super(Smooth, self).__init__()
        self.recursive = recursive

    def smooth(self, coords):
        '''
        Abstract method for smoothing method implementation.

        :param Iterable coords: Input coordinates to be smoothed.
        '''
        raise NotImplementedError("This is an abstract method.")

    def __call__(self, coords):
        '''
        Call method for all smoothing methods.

        Input coordinates should be iterable and each element should be numpy.ndarray. If length of :attr:`coords` is less then 3 smoothing method is not run and coordinates are returned unchanged.

        If :attr:`recursive` is set smoothing method is applied appropriate number of times.

        :param Iterable coords: Input coordinates to be smoothed.
        :rtype: numpy.ndarray
        :return: Smoothed coordinates.
        '''
        if len(coords) < 3:
            # this make no sense if coords length is less than 3
            return coords
        if self.recursive:
            coords_smooth = None
            for r in xrange(self.recursive):
                if coords_smooth is not None:
                    coords_smooth = self.smooth(coords_smooth)
                else:
                    coords_smooth = self.smooth(coords)
            return coords_smooth
        return self.smooth(coords)


class GeneralWindow(object):
    '''
    Base class for window based smoothing methods.
    '''
    def __init__(self, function=np.mean, **kwargs):
        '''
        :param function function: Function to be used for averaging coordinates within a window.
        '''
        super(GeneralWindow, self).__init__()
        self.function = function

    @staticmethod
    def max_window_at_pos(pos, size):
        '''
        Method returns maximal possible window at given position of the list with given size of the list. Returned window fits in to the list of given size and is symmetrical.

        :param int pos: Position in question.
        :param int size: Length of the list.
        :rtype: 2 element tuple of int
        :return: Lowest possible bound and highest possible bound of the window.
        '''
        # size is length
        # pos is zero based
        min_dist_to_edge = min((pos - 0, size - pos - 1))
        if (pos - 0) <= (size - pos - 1):
            # first half case
            lo = 0
            hi = pos + min_dist_to_edge + 1
        else:
            lo = pos - min_dist_to_edge
            hi = size
        return lo, hi

    def check_bounds_at_max_window_at_pos(self, lb, ub, pos, size):
        '''
        Method checks if window fits in to maximal possible window calculated according to :meth:`max_window_at_pos`. If not window is corrected.

        :param int lb: Lower bound of the window in question.
        :param int ub: Upper bound of the window in question.
        :param int pos: Position in question.
        :param int size: Length of the list.
        :rtype: 2 element tuple of int
        :return: Lowest possible bound and highest possible bound of the window corrected to maximal possible window.
        '''
        assert ub > pos
        lo, hi = self.max_window_at_pos(pos, size)
        if (pos - 0) <= (size - pos - 1):
            # first half case
            lo = lb
            hi = min(ub, hi)
        else:
            lo = max(lb, lo)
            hi = ub
        return lo, hi


class IntWindow(GeneralWindow):
    '''
    Base class for all window smoothing methods that require integer window.
    '''
    def __init__(self, window=5, **kwargs):
        '''
        :param int window: One side size of the window.
        '''
        super(IntWindow, self).__init__(**kwargs)
        self.window = int(window)


class FloatWindow(GeneralWindow):
    '''
    Base class for all window smoothing methods that require float window.
    '''
    def __init__(self, window=5., **kwargs):
        '''
        :param float window: Size of the window.
        '''
        super(FloatWindow, self).__init__(**kwargs)
        self.window = float(window)


class WindowSmooth(Smooth, IntWindow):
    '''
    Defined size window smoothing.

    For each coordinate a symmetrical (if possible) window of size defined by :attr:`window` is created.
    In case of coordinates at the edges created window is truncated to the edges. Next, all coordinates within the window are averaged with a function defined by :attr:`function`. Resulting value(s) are the smoothed coordinates.
    '''

    def __init__(self, **kwargs):
        super(WindowSmooth, self).__init__(**kwargs)

    @arrayify
    def smooth(self, coords):
        '''
        :param Iterable coords: Input coordinates to be smoothed.
        '''
        n = len(coords)
        for pos in xrange(n):
            lo, hi = self.max_window_at_pos(pos, n)
            lo = max(pos - self.window, lo)
            hi = min(pos + self.window, hi)
            lo, hi = self.check_bounds_at_max_window_at_pos(lo, hi, pos, n)

            yield self.function(coords[slice(lo, hi, None)], 0).tolist()


class DistanceWindowSmooth(Smooth, FloatWindow):
    '''
    Distance defined size window smoothing.

    This is modification of :class:`WindowSmooth` method.
    The difference is in the definition of the window size. Here, it is an average distance between points of input coordinates. Thus, before smoothing average distance between all points is calculated and this value is used to calculate actual window size.

    Next, for each coordinate a symmetrical (if possible) window of size calculated in the first step is created.
    In case of coordinates at the edges created window is truncated to the edges. Next, all coordinates within the window are averaged with a function defined by :attr:`function`. Resulting value(s) are the smoothed coordinates.
    '''

    def __init__(self, **kwargs):
        super(DistanceWindowSmooth, self).__init__(**kwargs)

    @arrayify
    def smooth(self, coords):
        '''
        :param Iterable coords: Input coordinates to be smoothed.
        '''
        n = len(coords)

        cdiff = traces.diff(coords)
        avgw = self.function(cdiff)
        if avgw != 0:
            window = int(np.ceil(self.window / avgw))
        else:
            window = len(coords) / 10  # FIXME: Magic number!
        if window < 1:
            window = 1

        for pos in xrange(n):
            lo, hi = self.max_window_at_pos(pos, n)
            lo = max(pos - window, lo)
            hi = min(pos + window, hi)
            lo, hi = self.check_bounds_at_max_window_at_pos(lo, hi, pos, n)

            yield self.function(coords[slice(lo, hi, None)], 0).tolist()


class ActiveWindowSmooth(Smooth, FloatWindow):
    '''
    Active size window smoothing.

    Similarly to :class:`DistanceWindowSmooth` method the window size is defined as a distance. The difference is that the actual window size is calculated for each point separately. Thus, for each coordinate the window is calculated by examining the distance differences between points. In this method window is not necessarily symmetrical. Once window is calculated all coordinates within the window are averaged with a function defined by :attr:`function`. Resulting value(s) are the smoothed coordinates.
    '''
    def __init__(self, **kwargs):
        super(ActiveWindowSmooth, self).__init__(**kwargs)

    @arrayify
    def smooth(self, coords):
        '''
        :param Iterable coords: Input coordinates to be smoothed.
        '''
        n = len(coords)
        d = traces.diff(coords)

        for pos in xrange(n):
            # lower boudary
            lb = pos
            ld = 0
            while (lb > 0) and (ld < self.window):
                ld += d[lb - 1]
                lb -= 1
            # upper boundary
            ub = pos
            ud = 0
            while (ub < n) and (ud < self.window):
                if ub < n - 1:
                    ud += d[ub]
                ub += 1

            lo, hi = self.check_bounds_at_max_window_at_pos(lb, ub, pos, n)

            yield self.function(coords[slice(lo, hi, None)], 0).tolist()


class MaxStepSmooth(Smooth):
    '''
    Maximal step smoothing.

    This method moves thorough coordinates and calculates distance over the traversed path. If it is then :attr:`step` the coordinate is used as a "cardinal point". The beginning and the end of the path are also added to the list of cardinal points. Next, all cardinal points and points of linear interpolation between cardinal points are returned as smoothed coordinates. Number of interpolated points is in accordance to points skipped between cardinal points.
    '''
    def __init__(self, step=1., **kwargs):
        super(MaxStepSmooth, self).__init__(**kwargs)
        self.step = step

    @arrayify
    def smooth(self, coords):
        '''
        :param Iterable coords: Input coordinates to be smoothed.
        '''
        n = len(coords)

        for pos in xrange(n):
            current_coord = coords[pos]
            if pos == 0:
                # yield first
                yield current_coord
                last_coord = current_coord
                to_yield_count = 0
            else:
                if pos == n - 1:
                    # yield last
                    if to_yield_count:
                        for coord in traces.tracepoints(last_coord, current_coord, to_yield_count):
                            yield coord
                    yield current_coord
                    last_coord = current_coord
                    to_yield_count = 0
                else:
                    current_step = traces.diff(np.vstack((current_coord, last_coord)))
                    if current_step > self.step:
                        # yield next!
                        if to_yield_count:
                            for coord in traces.tracepoints(last_coord, current_coord, to_yield_count):
                                yield coord
                        yield current_coord
                        last_coord = current_coord
                        to_yield_count = 0
                    else:
                        # update to yield count
                        to_yield_count += 1


class SavgolSmooth(Smooth):
    '''
    Savitzky-Golay based smoothing.

    Method uses 1D filter available in SciPy, see  :func:`~scipy.signal.savgol_filter`.
    For each dimension filter is applied separately. Only :attr:`window_length` and :attr:`polyorder` attributes are used.
    '''
    def __init__(self,window_length=5, polyorder=2, **kwargs):
        '''
        :param: int window_length: Size of the window, odd number.
        :param: int polyorder: Polynomial order.
        '''
        super(SavgolSmooth, self).__init__(**kwargs)
        self.savgol = partial(savgol_filter, window_length=window_length,polyorder=polyorder, axis=0, **kwargs)

    @arrayify
    def smooth(self, coords):
        '''
        :param Iterable coords: Input coordinates to be smoothed.
        '''
        return make_default_array(self.savgol(make_default_array(coords)))


class WindowOverMaxStepSmooth(Smooth):
    '''
    Window smoothing over maximal step smoothing.

    First, :class:`MaxStepSmooth` is applied, and then :class:`WindowSmooth`.
    '''
    def __init__(self, **kwargs):
        super(WindowOverMaxStepSmooth, self).__init__(**kwargs)
        self.window = WindowSmooth(**kwargs)
        self.mss = MaxStepSmooth(**kwargs)

    def smooth(self, coords):
        '''
        :param Iterable coords: Input coordinates to be smoothed.
        '''
        return self.window.smooth(self.mss.smooth(coords))


class ActiveWindowOverMaxStepSmooth(Smooth):
    '''
    Active window smoothing over maximal step smoothing.

    First, :class:`MaxStepSmooth` is applied, and then :class:`ActiveWindowSmooth`.
    '''
    def __init__(self, **kwargs):
        super(ActiveWindowOverMaxStepSmooth, self).__init__(**kwargs)
        self.window = ActiveWindowSmooth(**kwargs)
        self.mss = MaxStepSmooth(**kwargs)

    def smooth(self, coords):
        '''
        :param Iterable coords: Input coordinates to be smoothed.
        '''
        return self.window.smooth(self.mss.smooth(coords))


class DistanceWindowOverMaxStepSmooth(Smooth):
    '''
    Distance window smoothing over maximal step smoothing.

    First, :class:`MaxStepSmooth` is applied, and then :class:`DistanceWindowSmooth`.
    '''
    def __init__(self, **kwargs):
        super(DistanceWindowOverMaxStepSmooth, self).__init__(**kwargs)
        self.window = DistanceWindowSmooth(**kwargs)
        self.mss = MaxStepSmooth(**kwargs)

    def smooth(self, coords):
        '''
        :param Iterable coords: Input coordinates to be smoothed.
        '''
        return self.window.smooth(self.mss.smooth(coords))
