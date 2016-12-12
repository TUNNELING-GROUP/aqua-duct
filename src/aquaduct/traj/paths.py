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
from scipy.spatial.distance import cdist, pdist
from aquaduct.geom import traces
from aquaduct.traj.inlets import Inlet, InletTypeCodes
from aquaduct.utils.helpers import is_number, lind
from aquaduct.utils.helpers import tupleify, sortify, listify, arrayify1
from aquaduct.utils.maths import make_default_array


########################################################################################################################
# paths/list manipulations
# following part of code comes directly from the very initial tcl/vmd implementation
# all following functions should return list

def union(a, b):
    return [aa for aa in a if aa in b]


def glue(a, b):
    if a[-1] >= b[0]:
        g = list(set(a + b))
        g.sort()
        return g
    return []


@listify
def xor(a, b):
    ab = union(a, b)
    for e in glue(a, b):
        if e not in ab:
            yield e


def left(a, b):
    return union(a, xor(a, b))


def right(a, b):
    return union(b, xor(a, b))


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

        self.types = []
        self.frames = []
        self.coords = []

        # following is required to correct in and out paths that begin or end in scope and
        # begin or end at the very begining of MD or at very end of MD

        self.max_possible_frame = max_pf
        self.min_possible_frame = min_pf

    def add_coord(self, coord):
        # store coord as numpy float 32
        self.coords.append(make_default_array(coord))

    def add_object(self, frame):
        self.add_type(frame, self.object_name)

    def add_scope(self, frame):
        self.add_type(frame, self.scope_name)

    def add_type(self, frame, ftype):
        self.types.append(ftype)
        self.frames.append(frame)

    @property
    def max_frame(self):
        return max(self.frames)

    @property
    def min_frame(self):
        return min(self.frames)

    # xrange(self.max_frame,self.min_frame-1,-1)

    def get_paths_in(self):
        frames_range = xrange(self.max_frame, self.min_frame - 1, -1)
        for path in self.get_paths_for_frames_range(frames_range):
            path.sort()
            yield path

    def get_paths_out(self):
        frames_range = xrange(self.min_frame, self.max_frame + 1, 1)
        for path in self.get_paths_for_frames_range(frames_range):
            path.sort()
            yield path

    @sortify
    def get_paths_for_frames_range(self, frames_range):
        paths = []
        current_path = []
        # reverse
        for f in frames_range:
            # detect type
            if f in self.frames:
                t = self.types[self.frames.index(f)]
            else:
                t = self.out_name
            # use type to do something
            if t == self.out_name:
                if len(current_path) > 0:
                    yield current_path
                    current_path = []
            else:
                if len(current_path) > 0:
                    current_path.append(f)
                else:
                    if t == self.object_name:
                        current_path.append(f)
        if len(current_path) > 0:
            yield current_path

    def find_paths(self, fullonly=False):
        paths_out = list(self.get_paths_out())
        paths_in = list(self.get_paths_in())

        for path_in in paths_in:
            path_out = paths_out[0]
            path_glue = glue(path_in, path_out)
            if len(path_glue) > 0:
                path_out = paths_out.pop(0)
                path_core = union(path_in, path_out)
                path_in = left(path_in, path_out)
                path_out = right(path_core, path_out)
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

    def find_paths_coords(self, fullonly=False):
        for path in self.find_paths(fullonly=fullonly):
            yield self.get_single_path_coords(path)

    def find_paths_types(self, fullonly=False):
        for path in self.find_paths(fullonly=fullonly):
            yield self.get_single_path_types(path)

    def find_paths_coords_types(self, fullonly=False):
        for path in self.find_paths(fullonly=fullonly):
            yield path, self.get_single_path_coords(path), self.get_single_path_types(path)

    def get_single_path_coords(self, spath):
        # returns coordinates for single path
        # single path comprises of in,scope,out parts

        p_in, p_object, p_out = spath

        if len(self.frames) == len(self.coords):
            # not full trajectory
            in_ = []
            for f in p_in:
                in_.append(self.coords[self.frames.index(f)])
            object_ = []
            for f in p_object:
                object_.append(self.coords[self.frames.index(f)])
            out_ = []
            for f in p_out:
                out_.append(self.coords[self.frames.index(f)])
        else:
            # full trajectory
            in_ = []
            for f in p_in:
                in_.append(self.coords[f])
            object_ = []
            for f in p_object:
                object_.append(self.coords[f])
            out_ = []
            for f in p_out:
                out_.append(self.coords[f])

        return in_, object_, out_

    def get_single_path_types(self, spath):
        # returns typess for single path
        # single path comprises of in,scope,out parts

        # TODO: join it with get_single_path_coords

        p_in, p_object, p_out = spath

        if len(self.frames) == len(self.types):
            # not full trajectory
            in_ = []
            for f in p_in:
                in_.append(self.types[self.frames.index(f)])
            object_ = []
            for f in p_object:
                object_.append(self.types[self.frames.index(f)])
            out_ = []
            for f in p_out:
                out_.append(self.types[self.frames.index(f)])
        else:
            # full trajectory
            in_ = []
            for f in p_in:
                in_.append(self.types[f])
            object_ = []
            for f in p_object:
                object_.append(self.types[f])
            out_ = []
            for f in p_out:
                out_.append(self.types[f])

        return in_, object_, out_

    def barber_with_spheres(self, spheres):
        for center, radius in spheres:
            if len(self.coords) > 0:
                distances = cdist(np.matrix(center), self.coords, metric='euclidean').flatten()
                self.coords = lind(self.coords, np.argwhere(distances > radius).flatten().tolist())
                self.types = lind(self.types, np.argwhere(distances > radius).flatten().tolist())
                self.frames = lind(self.frames, np.argwhere(distances > radius).flatten().tolist())


# SinglePathID = namedtuple('SinglePathID', 'id nr')
class SinglePathID(object):
    def __init__(self, path_id=None, nr=None):
        self.id = path_id
        self.nr = nr

    def __str__(self):
        return '%d:%d' % (self.id, self.nr)


def yield_single_paths(gps, fullonly=False, progress=False):
    # iterates over gps - list of GenericPaths objects and transforms them in to SinglePath objects
    nr_dict = {}
    for nr, gp in enumerate(gps):
        path_id = gp.id
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
        self.path_in, self.path_object, self.path_out = paths
        self.coords_in, self.coords_object, self.coords_out = map(make_default_array, coords)
        self.types_in, self.types_object, self.types_out = types

        self.smooth_coords_in, self.smooth_coords_object, self.smooth_coords_out = None, None, None
        self.smooth_method = None

        # return np.vstack([c for c in self._coords if len(c) > 0])

    # ---------------------------------------------------------------------------------------------------------------- #
    # DEPRECATED
    @property
    def coords_first_in(self):
        if len(self.path_in) > 0:
            return self.coords_in[0]

    @property
    def coords_last_out(self):
        if len(self.path_out) > 0:
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
        return ([self.path_in_code] * len(self.path_in),
                [self.path_object_code] * len(self.path_object),
                [self.path_out_code] * len(self.path_out))

    @property
    def types_cont(self):
        return ([self.path_in_code] * len(self.path_in)) + ([self.path_object_code] * len(self.path_object)) + (
            [self.path_out_code] * len(self.path_out))

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
        return sum(map(len, self.paths))

    @property
    def begins(self):
        return self.paths_cont[0]

    @property
    def ends(self):
        return self.paths_cont[-1]

    @property
    def has_in(self):
        return len(self.path_in) > 0

    @property
    def has_object(self):
        return len(self.path_object) > 0

    @property
    def has_out(self):
        return len(self.path_out) > 0

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
        for path in self.paths:
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
        length = make_default_array(np.hstack((np.array([0]), np.cumsum(traces.diff(self.get_coords_cont(smooth=smooth))))))
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
