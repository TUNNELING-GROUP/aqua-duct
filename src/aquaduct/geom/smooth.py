# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2016  Tomasz Magdziarz, Alicja Płuciennik, Michał Stolarczyk <info@aquaduct.pl>
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

from aquaduct.geom import traces
from aquaduct.utils.helpers import arrayify


class Smooth(object):
    def __init__(self, recursive=None, **kwargs):

        self.recursive = recursive

    def smooth(self, coords):
        raise NotImplementedError("This is an abstract method.")

    def __call__(self, coords):
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
    def __init__(self,function=np.mean):
        self.function = function

    @staticmethod
    def max_window_at_pos(pos, size):
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
    def __init__(self,window=5,**kwargs):
        super(IntWindow,self).__init__(**kwargs)
        self.window = int(window)

class FloatWindow(GeneralWindow):
    def __init__(self,window=5.,**kwargs):
        super(FloatWindow,self).__init__(**kwargs)
        self.window = float(window)


class WindowSmooth(Smooth, IntWindow):
    def __init__(self, **kwargs):
        super(WindowSmooth,self).__init__(**kwargs)

    @arrayify
    def smooth(self, coords):
        n = len(coords)

        for pos in xrange(n):
            lo, hi = self.max_window_at_pos(pos, n)
            lo = max(pos - self.window, lo)
            hi = min(pos + self.window, hi)
            lo, hi = self.check_bounds_at_max_window_at_pos(lo, hi, pos, n)

            yield self.function(coords[slice(lo, hi, None)], 0).tolist()


class DistanceWindowSmooth(Smooth, FloatWindow):
    def __init__(self, **kwargs):
        super(DistanceWindowSmooth,self).__init__(**kwargs)

    @arrayify
    def smooth(self, coords):
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
    def __init__(self, **kwargs):
        super(ActiveWindowSmooth,self).__init__(**kwargs)

    @arrayify
    def smooth(self, coords):

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
    def __init__(self, step=1., **kwargs):
        super(MaxStepSmooth,self).__init__(**kwargs)
        self.step = step

    @arrayify
    def smooth(self, coords):

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


class WindowOverMaxStepSmooth(Smooth):
    def __init__(self, **kwargs):
        super(WindowOverMaxStepSmooth,self).__init__(**kwargs)

        self.window = WindowSmooth(**kwargs)
        self.mss = MaxStepSmooth(**kwargs)

    def smooth(self, coords):
        return self.window.smooth(self.mss.smooth(coords))


class ActiveWindowOverMaxStepSmooth(Smooth):
    def __init__(self, **kwargs):
        super(ActiveWindowOverMaxStepSmooth,self).__init__(**kwargs)

        self.window = ActiveWindowSmooth(**kwargs)
        self.mss = MaxStepSmooth(**kwargs)

    def smooth(self, coords):
        return self.window.smooth(self.mss.smooth(coords))


class DistanceWindowOverMaxStepSmooth(Smooth):
    def __init__(self, **kwargs):
        super(DistanceWindowOverMaxStepSmooth,self).__init__(**kwargs)

        self.window = DistanceWindowSmooth(**kwargs)
        self.mss = MaxStepSmooth(**kwargs)

    def smooth(self, coords):
        return self.window.smooth(self.mss.smooth(coords))
