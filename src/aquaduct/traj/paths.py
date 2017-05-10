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

import logging

logger = logging.getLogger(__name__)

from aquaduct.utils import clui
import numpy as np
from scipy.spatial.distance import cdist, pdist
from aquaduct.geom import traces
from aquaduct.traj.inlets import Inlet, InletTypeCodes
from aquaduct.utils.helpers import is_number, lind, list_blocks_to_slices, split_list
from aquaduct.utils.helpers import tupleify, sortify, listify, arrayify1
from aquaduct.utils.maths import make_default_array


########################################################################################################################
# paths/list manipulations
# following part of code comes directly from the very initial tcl/vmd implementation
# all following functions should return list
# functions union_smartr and xor_smartr were added later in python implementation to speed up calculations

def union_full(a, b):
    return [aa for aa in a if aa in b]


def union_smartr(a, b):
    b_ = SmartRange(b)
    return [aa for aa in a if b_.isin(aa)]


def union(a, b, smartr=True):
    if smartr:
        return union_smartr(a, b)
    return union_full(a, b)


def glue(a, b):
    if a[-1] >= b[0]:
        g = list(set(a + b))
        g.sort()
        return g
    return []


@listify
def xor_full(a, b):
    ab = union(a, b)
    for e in glue(a, b):
        if e not in ab:
            yield e


@listify
def xor_smartr(a, b):
    ab = SmartRange(union_smartr(a, b))
    for e in glue(a, b):
        if not ab.isin(e):
            yield e


def xor(a, b, smartr=True):
    if smartr:
        return xor_smartr(a, b)
    return xor_full(a, b)


def left(a, b, smartr=True):
    return union(a, xor(a, b, smartr=smartr), smartr=smartr)


def right(a, b, smartr=True):
    return union(b, xor(a, b, smartr=smartr), smartr=smartr)


########################################################################################################################

class SmartRangeFunction(object):
    def __init__(self, element, times):
        self.element = element
        self.times = times

    def __repr__(self):
        return "%s(%r,%d)" % (self.__class__.__name__, self.element, self.times)

    def get(self):
        raise NotImplementedError('This method should be implemented in a child class.')

    def rev(self):
        raise NotImplementedError('This method should be implemented in a child class.')

    def isin(self, element):
        raise NotImplementedError('This method should be implemented in a child class.')


class SmartRangeEqual(SmartRangeFunction):
    def get(self):
        return [self.element] * self.times

    def rev(self):
        return self

    def isin(self, element):
        return element == self.element


class SmartRangeIncrement(SmartRangeFunction):
    def get(self):
        return (self.element + i for i in xrange(self.times))

    def rev(self):
        return SmartRangeDecrement(self.element + self.times - 1, self.times)

    def isin(self, element):
        return (element >= self.element) and (element <= self.element + self.times - 1)


class SmartRangeDecrement(SmartRangeFunction):
    def get(self):
        return (self.element - i for i in xrange(self.times))

    def rev(self):
        return SmartRangeIncrement(self.element - self.times + 1, self.times)

    def isin(self, element):
        return (element <= self.element) and (element >= self.element - self.times + 1)


class SmartRange(object):
    def __init__(self, iterable=None):
        self.__elements = []
        self.__len = 0
        self.__min = None
        self.__max = None

        if iterable is not None:
            map(self.append, iterable)

    def last_element(self):
        if len(self.__elements) == 0:
            return None
        element = self.__elements[-1]
        if isinstance(element, SmartRangeFunction):
            return element.element
        return element

    def last_times(self):
        if len(self.__elements) == 0:
            return 0
        element = self.__elements[-1]
        if isinstance(element, SmartRangeFunction):
            return element.times
        return 1

    @property
    @listify
    def raw(self):
        for element in self.__elements:
            if not isinstance(element, SmartRangeFunction):
                yield SmartRangeEqual(element, 1)
            else:
                yield element

    def append(self, element):
        assert not isinstance(element, SmartRangeFunction)
        if len(self.__elements) == 0:
            self.__elements.append(element)
            self.__min = element
            self.__max = element
        else:
            if element == self.last_element():
                if isinstance(self.__elements[-1], SmartRangeEqual) or (
                        not isinstance(self.__elements[-1], SmartRangeFunction)):
                    self.__elements[-1] = SmartRangeEqual(element, self.last_times() + 1)
                else:
                    self.__elements.append(element)
            else:
                if not is_number(element):
                    self.__elements.append(element)
                else:
                    if element - self.last_times() == self.last_element():
                        if isinstance(self.__elements[-1], SmartRangeIncrement) or (
                                not isinstance(self.__elements[-1], SmartRangeFunction)):
                            self.__elements[-1] = SmartRangeIncrement(self.last_element(), self.last_times() + 1)
                        else:
                            self.__elements.append(element)
                    elif element + self.last_times() == self.last_element():
                        if isinstance(self.__elements[-1], SmartRangeDecrement) or (
                                not isinstance(self.__elements[-1], SmartRangeFunction)):
                            self.__elements[-1] = SmartRangeDecrement(self.last_element(), self.last_times() + 1)
                        else:
                            self.__elements.append(element)
                    else:
                        self.__elements.append(element)
            if element > self.__max:
                self.__max = element
            if element < self.__min:
                self.__min = element
        self.__len += 1

    def get(self):
        for element in self.__elements:
            if not isinstance(element, SmartRangeFunction):
                yield element
            else:
                for e in element.get():
                    yield e

    def rev(self):
        elements = []
        for e in self.__elements[::-1]:
            if isinstance(e, SmartRangeFunction):
                elements.append(e.rev())
            else:
                elements.append(e)
        self.__elements = elements

    def __len__(self):
        return self.__len

    def min(self):
        return self.__min

    def max(self):
        return self.__max

    def isin(self, element):
        for block in self.raw:
            if block.isin(element):
                return True
        return False


########################################################################################################################

class PathTypesCodes():
    path_in_code = 'i'
    path_object_code = 'c'
    path_out_code = 'o'


class GenericPathTypeCodes():
    object_name = 'c'
    scope_name = 's'
    out_name = 'n'


class GenericPaths(object, GenericPathTypeCodes):
    # object to store paths... is it required?

    def __init__(self, id_of_res, min_pf=None, max_pf=None):

        # id is any type of object; it is used as identifier

        self.id = id_of_res
        self.__types = SmartRange()
        self.__frames = SmartRange()
        self.coords = []

        # following is required to correct in and out paths that begin or end in scope and
        # begin or end at the very begining of MD or at very end of MD

        self.max_possible_frame = max_pf
        self.min_possible_frame = min_pf

    # info methods

    @property
    def types(self):
        return list(self.__types.get())

    @property
    def frames(self):
        return list(self.__frames.get())

    @property
    def max_frame(self):
        return self.__frames.max()

    @property
    def min_frame(self):
        return self.__frames.min()

    # add methods

    def add_coord(self, coord):
        # store coord as numpy float 32
        self.coords.append(make_default_array(coord))

    def add_object(self, frame):
        self.add_type(frame, self.object_name)

    def add_scope(self, frame):
        self.add_type(frame, self.scope_name)

    def add_type(self, frame, ftype):
        self.__types.append(ftype)
        self.__frames.append(frame)

    def _gpo(self):
        n = len(self.__frames)
        types = self.types
        begin = 0
        for block in self.__frames.raw:
            end = begin + block.times
            # get types of this block
            block_frames = list(block.get())
            block_types = types[begin:end]
            begin = end
            # now iterate over block_types in a search of out_name
            if self.out_name in block_types:
                while self.out_name in block_types:
                    to_yield = block_frames[:block_types.index(self.out_name)]
                    to_yield_types = block_types[:block_types.index(self.out_name)]
                    if len(to_yield) > 0:
                        while to_yield_types[0] != self.object_name:
                            to_yield.pop(0)
                            to_yield_types.pop(0)
                            if len(to_yield) == 0:
                                break
                        if len(to_yield) > 0:
                            yield to_yield
                    block_types = block_types[block_types.index(self.out_name):]
                    block_frames = block_frames[block_types.index(self.out_name):]
            if len(block_frames) > 0:
                while block_types[0] != self.object_name:
                    block_frames.pop(0)
                    block_types.pop(0)
                    if len(block_frames) == 0:
                        break
                if len(block_frames) > 0:
                    yield block_frames

    def _gpi(self):
        n = len(self.__frames)
        types = self.types
        begin = 0
        for block in self.__frames.raw:
            end = begin + block.times
            # get types of this block
            block_frames = list(block.get())
            block_types = types[begin:end]
            begin = end
            # now iterate over block_types in a search of out_name
            if self.out_name in block_types:
                while self.out_name in block_types:
                    to_yield = block_frames[:block_types.index(self.out_name)]
                    to_yield_types = block_types[:block_types.index(self.out_name)]
                    if len(to_yield) > 0:
                        while to_yield_types[-1] != self.object_name:
                            to_yield.pop(-1)
                            to_yield_types.pop(-1)
                            if len(to_yield) == 0:
                                break
                        if len(to_yield) > 0:
                            yield to_yield
                    block_types = block_types[block_types.index(self.out_name):]
                    block_frames = block_frames[block_types.index(self.out_name):]
            if len(block_frames) > 0:
                while block_types[-1] != self.object_name:
                    block_frames.pop(-1)
                    block_types.pop(-1)
                    if len(block_frames) == 0:
                        break
                if len(block_frames) > 0:
                    yield block_frames

    def get_paths_in(self):
        return self._gpi()

    def get_paths_out(self):
        return self._gpo()

    def find_paths(self, fullonly=False, smartr=True):
        paths_out = list(self.get_paths_out())
        paths_in = list(self.get_paths_in())

        for path_in in paths_in:
            path_out = paths_out[0]
            path_glue = glue(path_in, path_out)
            if len(path_glue) > 0:
                path_out = paths_out.pop(0)
                path_core = union(path_in, path_out, smartr=smartr)
                path_in = left(path_in, path_out, smartr=smartr)
                path_out = right(path_core, path_out, smartr=smartr)
            else:
                path_core = []
                path_out = []
            if len(path_in) == 0 or len(path_core) == 0 or len(path_out) == 0:
                if fullonly:
                    continue
            # TODO: make separate method/function for that
            # correct for max/min possible frame
            if self.max_possible_frame is not None:
                if len(path_out) > 0:
                    if path_out[-1] >= self.max_possible_frame:
                        path_core += path_out
                        path_out = []
            if self.min_possible_frame is not None:
                if len(path_in) > 0:
                    if path_in[0] <= self.min_possible_frame:
                        path_core = path_in + path_core
                        path_in = []
            yield path_in, path_core, path_out
        if not fullonly:
            path_in = []
            path_core = []
            for path_out in paths_out:
                # correct for max/min possible frame
                if self.max_possible_frame is not None:
                    if len(path_out) > 0:
                        if path_out[-1] >= self.max_possible_frame:
                            path_core += path_out
                            path_out = []
                yield path_in, path_core, path_out

    def find_paths_coords_types(self, fullonly=False):
        for path in self.find_paths(fullonly=fullonly):
            # yield path, self.get_single_path_coords(path), self.get_single_path_types(path)
            coords, types = self.get_single_path_coords_types(path)
            yield path, coords, types

    def get_single_path_coords_types(self, spath):
        # returns typess for single path
        # single path comprises of in,scope,out parts

        # TODO: join it with get_single_path_coords

        p_in, p_object, p_out = spath

        frames = self.frames
        types = self.types

        def get_be(p_):
            # if len(p_) == 1:
            # quit    return 0,1
            return frames.index(p_[0]), frames.index(p_[-1]) + 1

        in_t = []
        object_t = []
        out_t = []
        in_c = []
        object_c = []
        out_c = []

        if len(frames) == len(types):
            # not full trajectory
            # p_in
            if len(p_in) > 0:
                b, e = get_be(p_in)
                in_t = types[b:e]
                in_c = self.coords[b:e]
            # p_object
            if len(p_object) > 0:
                b, e = get_be(p_object)
                object_t = types[b:e]
                object_c = self.coords[b:e]
            # p_out
            if len(p_out) > 0:
                b, e = get_be(p_out)
                out_t = types[b:e]
                out_c = self.coords[b:e]
        else:
            # full trajectory
            # p_in
            if len(p_in) > 0:
                b, e = p_in[0], p_in[-1] + 1
                in_t = types[b:e]
                in_c = self.coords[b:e]
            # p_object
            if len(p_object) > 0:
                b, e = p_object[0], p_object[-1] + 1
                object_t = types[b:e]
                object_c = self.coords[b:e]
            # p_out
            if len(p_out) > 0:
                b, e = p_out[0], p_out[-1] + 1
                out_t = types[b:e]
                out_c = self.coords[b:e]

        return (in_c, object_c, out_c), (in_t, object_t, out_t)

    def barber_with_spheres(self, spheres):
        # calculate big distance matrix
        # distances = cdist(self.coords, [s.center for s in spheres],metric='euclidean')
        # compare with radii
        # tokeep = distances > np.matrix([[s.radius for s in spheres]])
        if len(spheres):
            tokeep = np.argwhere((cdist(self.coords, [s.center for s in spheres], metric='euclidean') > np.matrix(
                [[s.radius for s in spheres]])).all(1).A1).flatten().tolist()
            #tokeep = np.argwhere(tokeep).flatten().tolist()
            self.coords = lind(self.coords, tokeep)
            self.__types = SmartRange(lind(self.types, tokeep))
            self.__frames = SmartRange(lind(self.frames, tokeep))


# SinglePathID = namedtuple('SinglePathID', 'id nr')
class SinglePathID(object):
    def __init__(self, path_id=None, nr=None):
        self.id = path_id
        self.nr = nr

    def __str__(self):
        return '%d:%d' % (self.id, self.nr)

    def __eq__(self, other):
        return str(self) == str(other)

    def __ne__(self, other):
        return str(self) != str(other)

def yield_single_paths(gps, fullonly=False, progress=False):
    # iterates over gps - list of GenericPaths objects and transforms them in to SinglePath objects
    nr_dict = {}
    for nr, gp in enumerate(gps):
        path_id = gp.id
        with clui.tictoc('Processing path %d' % path_id):
            for paths, coords, types in gp.find_paths_coords_types(fullonly=fullonly):
                if path_id in nr_dict:
                    nr_dict.update({path_id: (nr_dict[path_id] + 1)})
                else:
                    nr_dict.update({path_id: 0})
                if progress:
                    yield SinglePath(SinglePathID(path_id=path_id, nr=nr_dict[path_id]), paths, coords, types), nr
                else:
                    yield SinglePath(SinglePathID(path_id=path_id, nr=nr_dict[path_id]), paths, coords, types)


class SinglePath(object, PathTypesCodes, InletTypeCodes):
    # special class
    # represents one path

    empty_coords = make_default_array(np.zeros((0, 3)))

    def __init__(self, path_id, paths, coords, types):

        self.id = path_id
        # for paths use SmartRanges and then provide methods to read them
        self.__path_in, self.__path_object, self.__path_out = map(SmartRange, paths)
        # similarly, do it with types
        # self.path_in, self.path_object, self.path_out = paths
        self.__types_in, self.__types_object, self.__types_out = map(SmartRange, types)
        # self.types_in, self.types_object, self.types_out = types

        self.coords_in, self.coords_object, self.coords_out = map(make_default_array, coords)

        self.smooth_coords_in, self.smooth_coords_object, self.smooth_coords_out = None, None, None
        self.smooth_method = None

        # return np.vstack([c for c in self._coords if len(c) > 0])

    @property
    def path_in(self):
        return list(self.__path_in.get())

    @property
    def path_object(self):
        return list(self.__path_object.get())

    @property
    def path_out(self):
        return list(self.__path_out.get())

    @property
    def types_in(self):
        return list(self.__types_in.get())

    @property
    def types_object(self):
        return list(self.__types_object.get())

    @property
    def types_out(self):
        return list(self.__types_out.get())

    # ---------------------------------------------------------------------------------------------------------------- #
    # DEPRECATED
    @property
    def coords_first_in(self):
        if len(self.__path_in) > 0:
            return self.coords_in[0]

    @property
    def coords_last_out(self):
        if len(self.__path_out) > 0:
            return self.coords_out[-1]

    @property
    def coords_filo(self):
        # first in and last out plus type!
        return [(inlet, {0: self.incoming, 1: self.outgoing}[nr]) for nr, inlet in
                enumerate((self.coords_first_in, self.coords_last_out)) if inlet is not None]

    # ---------------------------------------------------------------------------------------------------------------- #

    def get_inlets(self):
        if self.has_in:
            yield Inlet(coords=self.coords_in[0],
                        type=(InletTypeCodes.surface, InletTypeCodes.incoming),
                        reference=self.id)
            yield Inlet(coords=self.coords_in[-1],
                        type=(InletTypeCodes.internal, InletTypeCodes.incoming),
                        reference=self.id)
        if self.has_out:
            yield Inlet(coords=self.coords_out[0],
                        type=(InletTypeCodes.internal, InletTypeCodes.outgoing),
                        reference=self.id)
            yield Inlet(coords=self.coords_out[-1],
                        type=(InletTypeCodes.surface, InletTypeCodes.outgoing),
                        reference=self.id)

    ####################################################################################################################
    # coords

    @property
    def coords(self):
        return self.coords_in, self.coords_object, self.coords_out

    @property
    def coords_cont(self):
        # returns coords as one array
        return make_default_array(np.vstack([c for c in self.coords if len(c) > 0]))

    ####################################################################################################################
    # paths

    @property
    def __paths(self):
        return self.__path_in, self.__path_object, self.__path_out

    @property
    def paths(self):
        return self.path_in, self.path_object, self.path_out

    @property
    def paths_cont(self):
        return self.path_in + self.path_object + self.path_out

    ####################################################################################################################
    # types

    @property
    def types(self):
        # spath types
        return ([self.path_in_code] * len(self.__path_in),
                [self.path_object_code] * len(self.__path_object),
                [self.path_out_code] * len(self.__path_out))

    @property
    def types_cont(self):
        return ([self.path_in_code] * len(self.__path_in)) + ([self.path_object_code] * len(self.__path_object)) + (
            [self.path_out_code] * len(self.__path_out))

    @property
    def gtypes(self):
        # generic types
        return self.types_in, self.types_object, self.types_out

    @property
    def gtypes_cont(self):
        return self.types_in + self.types_object + self.types_out

    @property
    @tupleify
    def etypes(self):
        # extended types
        for t, g in zip(self.types, self.gtypes):
            yield [''.join(t) for t in zip(t, g)]

    @property
    def etypes_cont(self):
        return [''.join(t) for t in zip(self.types_cont, self.gtypes_cont)]

    ####################################################################################################################

    @property
    def size(self):
        return sum(map(len, self.__paths))

    @property
    def begins(self):
        return self.paths_cont[0]

    @property
    def ends(self):
        return self.paths_cont[-1]

    @property
    def has_in(self):
        return len(self.__path_in) > 0

    @property
    def has_object(self):
        return len(self.__path_object) > 0

    @property
    def has_out(self):
        return len(self.__path_out) > 0

    ####################################################################################################################

    @tupleify
    def get_coords(self, smooth=None):
        # TODO: it is not used to get smooth coords but to get coords in general, conditionally smoothed
        # if smooth is not none applies smoothing
        if smooth is not None:
            if smooth != self.smooth_method:
                self.smooth_coords_in, self.smooth_coords_object, self.smooth_coords_out = self._make_smooth_coords(
                    self.coords_cont, smooth)
                self.smooth_method = smooth
            for nr, coords in enumerate((self.smooth_coords_in, self.smooth_coords_object, self.smooth_coords_out)):
                if coords is None:
                    yield self.coords[nr]
                else:
                    yield coords
        else:
            for coords in self.coords:
                yield coords

    def get_coords_cont(self, smooth=None):
        # returns coords as one array
        return make_default_array(np.vstack([c for c in self.get_coords(smooth) if len(c) > 0]))

    @tupleify
    def _make_smooth_coords(self, coords, smooth=None):
        # if smooth is not none applies smoothing
        # smooth should be callable and should return an object of length equal to submitted one
        # get continuous coords
        if smooth:
            coords_smooth = smooth(coords)
        else:
            coords_smooth = self.coords
        # now lets return tupple of coords
        nr = 0
        for path in self.__paths:
            if len(path) > 0:
                yield make_default_array(coords_smooth[nr:nr + len(path)])
                nr += len(path)
            else:
                yield self.empty_coords

    def apply_smoothing(self, smooth):
        # permament change!
        self.coords_in, self.coords_object, self.coords_out = self._make_smooth_coords(self.coords_cont, smooth)

    ####################################################################################################################

    def get_distance_cont(self, smooth=None, normalize=False):
        length = make_default_array(
            np.hstack((np.array([0]), np.cumsum(traces.diff(self.get_coords_cont(smooth=smooth))))))
        if normalize:
            if is_number(normalize):
                norm_factor = float(normalize)
            else:
                norm_factor = max(length)
            length = length / norm_factor
        return length

    def get_distance_rev_cont(self, *args, **kwargs):
        length = self.get_distance_cont(*args, **kwargs)
        return length[-1] - length

    @arrayify1
    def get_distance_both_cont(self, *args, **kwargs):
        length = self.get_distance_cont(*args, **kwargs)
        return (min(n, r) for n, r in zip(length, length[-1] - length))

    @arrayify1
    def get_velocity_cont(self, smooth=None, normalize=False):
        distance = self.get_distance_cont(smooth=smooth, normalize=normalize)
        return traces.derrivative(distance)

    @arrayify1
    def get_acceleration_cont(self, smooth=None, normalize=False):
        velocity = self.get_velocity_cont(smooth=smooth, normalize=normalize)
        return traces.derrivative(velocity)

        ####################################################################################################################


class MasterPath(SinglePath):
    def __init__(self, sp):
        SinglePath.__init__(self, sp.id, sp.paths, sp.coords, sp.gtypes)
        self.width_cont = None

    def add_width(self, width):
        assert len(width) == self.size
        self.width_cont = width
