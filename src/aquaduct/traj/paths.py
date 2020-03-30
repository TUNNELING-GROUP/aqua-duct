# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2016-2018  Tomasz Magdziarz, Alicja Płuciennik, Michał Stolarczyk <info@aquaduct.pl>
# Copyright (C) 2019  Tomasz Magdziarz <info@aquaduct.pl>
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


from collections import OrderedDict
from aquaduct.utils import clui
import numpy as np
from scipy.spatial.distance import cdist
from aquaduct.geom import traces
from aquaduct.traj.inlets import Inlet, InletTypeCodes
from aquaduct.utils.helpers import is_number, lind, glind, \
    SmartRange  # smart ranges are required here to provide backward compatibility with v0.3
from aquaduct.utils.sets import intersection, glue, left, right
from aquaduct.utils.helpers import tupleify, listify, arrayify1, list_blocks_to_slices
from aquaduct.utils.maths import make_default_array
from aquaduct.traj.sandwich import Reader, SingleResidueSelection



from array import array

from itertools import chain


class PathTypesCodes(object):
    __slots__ = ()  # "path_in_code path_object_code path_out_code path_walk_code".split()
    # path_in_code = 'i'
    # path_object_code = 'c'
    # path_out_code = 'o'
    # path_walk_code = 'w'
    path_in_code = 1
    path_object_code = 2
    path_out_code = 3
    path_walk_code = 0


class GenericPathTypeCodes(object):
    __slots__ = ()  # "object_name scope_name out_name".split()
    # object_name = 'c'
    # scope_name = 's'
    # out_name = 'n'
    object_name = 8
    scope_name = 9
    out_name = 10


class GenericPaths(GenericPathTypeCodes):
    # object to store paths... is it required?
    __slots__ = 'id single_res_selection name _types _frames max_possible_frame min_possible_frame'.split()

    def __init__(self, id_of_res, name_of_res=None,
                 min_pf=None, max_pf=None):

        super(GenericPaths, self).__init__()

        # id is any type of object; it is used as identifier
        # single_res_selection is object which have coords method that accepts frames and returns coordinates

        self.id = id_of_res
        self.single_res_selection = SingleResidueSelection(self.id)

        if name_of_res is not None:
            self.name = name_of_res
        else:
            self.name = 'UNK'  # FIXME: magic constant
        assert isinstance(self.name, str)
        assert len(self.name) == 3
        self._types = SmartRange()
        # self._frames = SmartRange()
        self._frames = array('i')

        # following is required to correct in and out paths that begin or end in scope and
        # begin or end at the very begining of MD or at very end of MD

        self.max_possible_frame = max_pf
        self.min_possible_frame = min_pf

    def update_types_frames(self, types, frames):
        if isinstance(types, SmartRange) and isinstance(frames, SmartRange):
            self._types = types
            self._frames = array('i', list(frames.get()))
        elif isinstance(types, SmartRange) and isinstance(frames, array):
            self._types = types
            self._frames = frames
        # make it from "any" iterable
        else:
            self._types = SmartRange(types)
            self._frames = array('i', frames)
        assert len(self._frames) == len(set(self._frames)), (len(self._frames), len(set(self._frames)))

    '''

    def __getstate__(self):
        return self.id, self.name, self._types, self._frames, self.max_possible_frame, self.min_possible_frame

    def __setstate__(self, state):
        # FIXME: tmp solution
        if isinstance(state, dict):
            self.id = state['id']
            self.name = state['name']
            self._types = state['_GenericPaths_types']
            self._frames = state['_GenericPaths_frames']
            self.max_possible_frame = state['max_possible_frame']
            self.min_possible_frame = state['min_possible_frame']
        else:
            self.id, self.name, self._types, self._frames, self.max_possible_frame, self.min_possible_frame = state

            Reader.reset()
        self.single_res_selection = SingleResidueSelection(self.id)

    '''

    # info methods
    @property
    def types(self):
        return list(self._types.get())

    @property
    def types_promise(self):
        return self._types.get()

    @property
    def frames_of_object(self):
        def get_foo():
            frame = self._frames.first_element()
            for sr in self._types.raw:
                if sr.element == self.object_name:
                    yield range(frame, frame + sr.times)
                frame += sr.times

        return chain(*get_foo())

    @property
    def frames_of_scope(self):
        def get_fos():
            frame = self._frames.first_element()
            for sr in self._types.raw:
                if sr.element == self.scope_name:
                    yield range(frame, frame + sr.times)
                frame += sr.times

        return chain(*get_fos())

    @property
    def frames(self):
        # return list(self._frames.get())
        return list(self._frames)

    @property
    def _frames_sr(self):
        return SmartRange(fast_array=self._frames)

    # @property
    # def frames_promise(self):
    #    return self._frames.get()

    def discard_singletons(self, singl=1, skiptype=GenericPathTypeCodes.object_name):
        # singl is chunk size to discard
        types = self.types
        new_types = []
        new_frames = []
        seek = 0
        # loop over frames chunks:
        for chunksr in self._frames_sr.raw:
            if len(chunksr) <= singl:
                if skiptype not in types[seek:seek + len(chunksr)]:
                    continue

            new_frames.append(chunksr)
            new_types.extend(types[seek:seek + len(chunksr)])
            seek += len(chunksr)

        self._types = SmartRange(new_types)
        self._frames = array('i', chain(*(ch.get() for ch in new_frames)))

    @property
    def coords(self):
        return self.single_res_selection.coords(self._frames_sr)

    @property
    def max_frame(self):
        return max(self._frames)

    @property
    def min_frame(self):
        return min(self._frames)

    # add methods
    def add_foos(self, foo, fos):
        foo_ = None
        fos_ = None
        while len(foo) + len(fos):
            if len(foo) and foo_ is None:
                foo_ = foo.pop(0)
            if len(fos) and fos_ is None:
                fos_ = fos.pop(0)
            if foo_ is None:
                self.add_scope(fos_)
                fos_ = None
            elif fos_ is None:
                self.add_object(foo_)
                foo_ = None
            elif foo_ > fos_:
                self.add_scope(fos_)
                fos_ = None
            elif fos_ > foo_:
                self.add_object(foo_)
                foo_ = None

    def add_012(self, os_in_frames, reset=False):
        if reset:
            self._types = SmartRange()
            self._frames = array('i')
        for frame, os_type in enumerate(os_in_frames):
            if os_type == 2:
                self.add_object(frame)
            elif os_type == 1:
                self.add_scope(frame)

    def add_object(self, frame):
        self.add_type(frame, self.object_name)

    def add_scope(self, frame):
        self.add_type(frame, self.scope_name)

    def add_type(self, frame, ftype):
        assert frame not in self._frames
        self._types.append(ftype)
        self._frames.append(frame)

    def add_frames_types(self, frames, types):
        for f, t in zip(frames, types):
            self._types.append(t)
            self._frames.append(f)

    '''
    def _split_path_by_edges_(self,path):
        # input: path tuple
        # output: path_cont
        edges = self.single_res_selection.get_edges()
        path_cont = []
        map(path_cont.extend, path)
        if edges:
            for e in edges:
                if e in path_cont:
                    path_cont_2y = path_cont[:path_cont.index(e) + 1]
                    path_cont = path_cont[path_cont.index(e) + 1:]
                    yield path_cont_2y
        if path_cont:
            yield path_cont
    '''

    '''
    def _split_path_by_types_(self,path):
        # input: path_cont
        # output: path tuple
        types = self.get_path_cont_types(path)
        # now, split it into single path or passing
        if self.object_name in types:
            # in
            bp = types.index(self.object_name)
            in_p = path[:bp]
            path = path[bp:]
            types = types[bp:]
            # obj
            types = types[::-1]
            path = path[::-1]
            ep = types.index(self.object_name)
            out_p = path[:ep]
            path = path[ep:]
            out_p = out_p[::-1]
            path = path[::-1]
            return in_p,path,out_p
        else:
            return (path,[],[]) # passing
    '''

    '''
    def _consider_edges(self,path,passing=None):
        for path_cont in self._split_path_by_edges_(path):
            spath = self._split_path_by_types_(path_cont)
            if len(spath[1])>0 or passing:
                yield spath
    '''

    def _gpt(self):
        # get, I'm just passing through
        n = len(self._frames)
        types = self.types
        begin = 0
        for block in self._frames_sr.raw:
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
                        if self.object_name not in to_yield_types:
                            yield to_yield
                    block_types = block_types[block_types.index(self.out_name):]
                    block_frames = block_frames[block_types.index(self.out_name):]
            if len(block_frames) > 0:
                if len(block_frames) > 0:
                    if self.object_name not in block_types:
                        yield block_frames

    def _gpo(self, frames_sr):
        n = len(frames_sr)
        types = self.types
        begin = 0
        for block in frames_sr.raw:
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

    def _gpi(self, frames_sr):
        n = len(frames_sr)
        types = self.types
        begin = 0
        for block in frames_sr.raw:
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

    def find_paths(self, fullonly=False):
        # this looks for normal, ie containing 'core' paths
        frames_sr = self._frames_sr
        paths_out = list(self._gpo(frames_sr))
        paths_in = list(self._gpi(frames_sr))
        del frames_sr

        for path_in in paths_in:
            path_out = paths_out[0]
            path_glue = glue(path_in, path_out)
            if len(path_glue) > 0:
                path_out = paths_out.pop(0)
                path_core = intersection(path_in, path_out)
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

    def find_paths_types(self, fullonly=False, passing=None):
        # both loops below yields correct paths
        # if waterfall let's split each of yielded paths again
        edges = self.single_res_selection.get_edges()
        for path in self.find_paths(fullonly=fullonly):
            '''
            for path_sbe in self._consider_edges(path,passing):
                types = self.get_single_path_types(path_sbe)
                # yield path, coords, types
                yield path_sbe, types
            '''
            types = self.get_single_path_types(path)
            # yield path, coords, types
            yield path, types
        if passing:
            for path in self._gpt():
                '''
                for path_sbe in self._consider_edges((path, [], []), passing):
                    types = self.get_single_path_types(path_sbe)
                    # yield path, coords, types
                    yield path_sbe, types
                '''
                types = self.get_single_path_types((path, [], []))
                # yield path, coords, types
                yield (path, [], []), types

    def get_path_cont_types(self, path_cont):

        frames = self.frames
        types = self.types

        def get_be(p_):
            # if len(p_) == 1:
            # quit    return 0,1
            return frames.index(p_[0]), frames.index(p_[-1]) + 1

        if len(frames) == len(types):
            # not full trajectory
            if len(path_cont) > 0:
                b, e = get_be(path_cont)
                return types[b:e]
        else:
            # full trajectory
            if len(path_cont) > 0:
                b, e = path_cont[0], path_cont[-1] + 1
                return types[b:e]
        return []

    def get_single_path_types(self, spath):
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
                # in_c = self.coords[b:e]
            # p_object
            if len(p_object) > 0:
                b, e = get_be(p_object)
                object_t = types[b:e]
                # object_c = self.coords[b:e]
            # p_out
            if len(p_out) > 0:
                b, e = get_be(p_out)
                out_t = types[b:e]
                # out_c = self.coords[b:e]
        else:
            # full trajectory
            # p_in
            if len(p_in) > 0:
                b, e = p_in[0], p_in[-1] + 1
                in_t = types[b:e]
                # in_c = self.coords[b:e]
            # p_object
            if len(p_object) > 0:
                b, e = p_object[0], p_object[-1] + 1
                object_t = types[b:e]
                # object_c = self.coords[b:e]
            # p_out
            if len(p_out) > 0:
                b, e = p_out[0], p_out[-1] + 1
                out_t = types[b:e]
                # out_c = self.coords[b:e]

        # return (in_c, object_c, out_c), (in_t, object_t, out_t)
        return in_t, object_t, out_t

    def barber_with_spheres(self, spheres):
        # calculate big distance matrix
        # distances = cdist(self.coords, [s.center for s in spheres],metric='euclidean')
        # compare with radii
        # tokeep = distances > np.matrix([[s.radius for s in spheres]])
        if len(spheres):
            tokeep = np.argwhere(
                (cdist(np.array(list(self.coords)), [s.center for s in spheres], metric='euclidean') > np.matrix(
                    [[s.radius for s in spheres]])).all(1).A1).flatten().tolist()
            # tokeep = np.argwhere(tokeep).flatten().tolist()
            self._types = SmartRange(lind(self.types, tokeep))
            self._frames = SmartRange(lind(self.frames, tokeep))


# SinglePathID = namedtuple('SinglePathID', 'id nr')
class SinglePathID(object):
    __slots__ = 'id nr name'.split()

    def __init__(self, path_id=None, nr=None, name=None):
        assert path_id is not None, "path_id cannot be None."
        self.id = path_id
        assert nr is not None, "nr cannot be None."
        self.nr = nr
        assert name is not None, "name cannot be None."
        self.name = name

    def __getstate__(self):
        return self.id, self.nr, self.name

    def __setstate__(self, state):
        self.id, self.nr, self.name = state

    def __str__(self):
        # by default name is not returned
        return '%d:%d:%d' % (self.id + (self.nr,))

    def __eq__(self, other):
        if isinstance(other, SinglePathID):
            return self.id == other.id and self.nr == other.nr and self.name == other.name
        return False

    def __hash__(self):
        return hash(str(self))

@listify
def yield_single_paths(gps, fullonly=None, progress=None, passing=None):
    # iterates over gps - list of GenericPaths objects and transforms them in to SinglePath objects
    nr_dict = {}
    for nr, gp in enumerate(gps):
        path_id = gp.id
        path_name = gp.name
        with clui.tictoc('Processing path %d:%d' % path_id):
            for paths, types in gp.find_paths_types(fullonly=fullonly, passing=passing):
                # for paths, coords, types in gp.find_paths_types(fullonly=fullonly):
                if path_id in nr_dict:
                    nr_dict.update({path_id: (nr_dict[path_id] + 1)})
                else:
                    nr_dict.update({path_id: 0})

                pid = SinglePathID(path_id=path_id, nr=nr_dict[path_id], name=path_name)

                if len(paths[1]) > 0:
                    # if everything is OK
                    sp = SinglePath(pid, paths, types)
                else:
                    # this is passing through
                    if passing:
                        # import ipdb as pdb; pdb.set_trace()
                        sp = PassingPath(pid, paths, types)
                        sp.has_in = sp.paths_first_in > gp.min_possible_frame
                        sp.has_out = sp.paths_last_out < gp.max_possible_frame
                    else:
                        continue
                if progress:
                    yield sp, nr
                else:
                    yield sp


def correct_spaths_ids(spaths, pbar):
    seen = {}
    for sp in spaths:
        if sp.id.id not in seen:
            seen.update({sp.id.id: 0})
        else:
            seen[sp.id.id] += 1
        sp.id.nr = seen[sp.id.id]
        pbar.next()
        # print str(sp.id)


def yield_generic_paths(spaths, progress=None):
    rid_seen = OrderedDict()


    number_of_frames = Reader.number_of_frames(onelayer=True) - 1
    t_object = GenericPathTypeCodes.object_name
    for sp in spaths:
        current_rid = sp.id.id
        if current_rid not in rid_seen:
            rid_seen.update(
                {current_rid: GenericPaths(current_rid, name_of_res=sp.id.name, min_pf=0, max_pf=number_of_frames)})
        # TODO: following loop is not an optimal solution, it is better to add types and frames in one call
        for t, f in zip(sp.gtypes_cont, sp.paths_cont):
            if t == t_object:
                rid_seen[current_rid].add_object(f)
            else:
                rid_seen[current_rid].add_scope(f)
        if progress:
            progress.next()
    # because paths stores now frames as array and produces smartranges on demand with fast_array option
    # it is required to keep frames (and types) in order, otherwise smartranges are wrong
    for p in rid_seen.values():
        new_order = np.argsort(p.frames)
        p.update_types_frames(glind(p.types, new_order), glind(p.frames, new_order))
        progress.next()

    return list(rid_seen.values())


class MacroMolPath(PathTypesCodes, InletTypeCodes):
    # special class
    # represents one path

    __slots__ = "id _path_in _path_object _path_out _types_in _types_object _types_out _object_len single_res_selection".split()

    empty_coords = make_default_array(np.zeros((0, 3)))

    def __init__(self, path_id, paths, types):

        super(MacroMolPath, self).__init__()

        self.id = path_id
        self.single_res_selection = SingleResidueSelection(self.id.id)
        # for paths use SmartRanges and then provide methods to read them
        self._path_in, self._path_object, self._path_out = list(map(SmartRange, paths))
        # similarly, do it with types
        # self.path_in, self.path_object, self.path_out = paths
        self._types_in, self._types_object, self._types_out = list(map(SmartRange, types))
        # self.types_in, self.types_object, self.types_out = types

        # self.coords_in, self.coords_object, self.coords_out = map(make_default_array, coords)

        # self.smooth_coords_in, self.smooth_coords_object, self.smooth_coords_out = None, None, None
        # self.smooth_method = None

        # return np.vstack([c for c in self._coords if len(c) > 0])

        self._object_len = None


    
    def __getstate__(self):
        return self.id, self._path_in, self._path_object, self._path_out, self._types_in, self._types_object, self._types_out, self._object_len

    def __setstate__(self, state):
        self.id, self._path_in, self._path_object, self._path_out, self._types_in, self._types_object, self._types_out, self._object_len = state
        self.single_res_selection = SingleResidueSelection(self.id.id)

    

    def add_paths4(self, path_in, path_object, path_object_strict, path_out):
        # init empty path
        self.__path_in, self.__path_object, self.__path_out = list(map(SmartRange, (path_in, path_object, path_out)))
        self.__types_in = SmartRange([GenericPathTypeCodes.scope_name] * len(path_in))
        self.__types_out = SmartRange([GenericPathTypeCodes.scope_name] * len(path_out))

        tobj = np.array([GenericPathTypeCodes.scope_name] * len(path_object))
        tobj[[nr for nr, f in enumerate(path_object) if f in path_object_strict]] = GenericPathTypeCodes.object_name
        self.__types_object = SmartRange(tobj.tolist())

    def _object_len_calculate(self):
        for nr, real_coords in enumerate(traces.midpoints(self.coords)):
            if nr != 1:
                continue
            if len(real_coords) <= 1:
                return 0.
            return float(sum(traces.diff(real_coords)))

    @property
    def object_len(self):
        if self._object_len is None:
            self._object_len = self._object_len_calculate()
        return self._object_len

    def is_single(self):
        raise NotImplementedError("Implementation missing.")

    def is_passing(self):
        raise NotImplementedError("Implementation missing.")

    def is_frame_in(self, frame):
        return self._path_in.isin(frame)

    def is_frame_object(self, frame):
        return self._path_object.isin(frame)

    def is_frame_out(self, frame):
        return self._path_out.isin(frame)

    def is_frame_walk(self, frame):
        return self.is_frame_in(frame) or self.is_frame_object(frame) or self.is_frame_out(frame)

    @property
    def path_in(self):
        return list(self._path_in.get())

    @property
    def path_object(self):
        return list(self._path_object.get())

    @property
    def path_object_strict_len(self):
        # number of frames strictly in object
        return self.etypes[1].count(self.etypes[1][0])

    @property
    def path_out(self):
        return list(self._path_out.get())

    @property
    def types_in(self):
        return list(self._types_in.get())

    @property
    def types_object(self):
        return list(self._types_object.get())

    @property
    def types_out(self):
        return list(self._types_out.get())

    # ---------------------------------------------------------------------------------------------------------------- #

    @property
    def coords_first_in(self):
        if len(self._path_in) > 0:
            return self.coords_in[0]
            # return self.single_res_selection.coords([self._path_in.first_element()])

    @property
    def paths_first_in(self):
        if len(self._path_in) > 0:
            return self.path_in[0]

    @property
    def coords_last_out(self):
        if len(self._path_out) > 0:
            return self.coords_out[-1]
            # return self.single_res_selection.coords([self._path_out.last_element()])

    @property
    def paths_last_out(self):
        if len(self._path_out) > 0:
            return self.path_out[-1]

    @property
    def coords_filo(self):
        # first in and last out plus type!
        return [(inlet, {0: self.incoming, 1: self.outgoing}[nr]) for nr, inlet in
                enumerate((self.coords_first_in, self.coords_last_out)) if inlet is not None]

    # ---------------------------------------------------------------------------------------------------------------- #

    def get_inlets(self):
        if self.has_in:
            yield Inlet(coords=self.coords_in[0],
                        inlet_type=(InletTypeCodes.surface, InletTypeCodes.incoming),
                        reference=self.id,
                        frame=self.path_in[0])
            yield Inlet(coords=self.coords_in[-1],
                        inlet_type=(InletTypeCodes.internal, InletTypeCodes.incoming),
                        reference=self.id,
                        frame=self.path_in[-1])
        if self.has_out:
            yield Inlet(coords=self.coords_out[0],
                        inlet_type=(InletTypeCodes.internal, InletTypeCodes.outgoing),
                        reference=self.id,
                        frame=self.path_out[0])
            yield Inlet(coords=self.coords_out[-1],
                        inlet_type=(InletTypeCodes.surface, InletTypeCodes.outgoing),
                        reference=self.id,
                        frame=self.path_out[-1])

    def remove_inlet(self, inlet_type):
        # only surface type can be removed
        if InletTypeCodes.surface == inlet_type[0]:
            if self.has_in and InletTypeCodes.incoming == inlet_type[1]:
                self._path_object = SmartRange(self.path_in + self.path_object)
                self._path_in = SmartRange([])
                self._types_object = SmartRange(self.types_in + self.types_object)
                self._types_in = SmartRange([])
            elif self.has_out and InletTypeCodes.outgoing == inlet_type[1]:
                self._path_object = SmartRange(self.path_object + self.path_out)
                self._path_out = SmartRange([])
                self._types_object = SmartRange(self.types_object + self.types_out)
                self._types_out = SmartRange([])
        # TODO: raise exceptions or at least warnings if wrong type of inlet is used

    ####################################################################################################################
    # coords

    @property
    def coords_in(self):
        return self.single_res_selection.coords(self._path_in)

    @property
    def coords_object(self):
        return self.single_res_selection.coords(self._path_object)

    @property
    def coords_object_strict(self):
        i = (nr for nr, et in enumerate(self.etypes[1]) if et == self.etypes[1][0])
        return self.coords_object[list(i)]

    @property
    def center_of_object(self):
        return np.mean(self.coords_object_strict, 0)

    @property
    def coords_out(self):
        return self.single_res_selection.coords(self._path_out)

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
    def _paths(self):
        return self._path_in, self._path_object, self._path_out

    @property
    def paths(self):
        return self.path_in, self.path_object, self.path_out

    @property
    def paths_cont(self):
        pathsc = []
        for p in self.paths:
            pathsc += p
        return pathsc

    ####################################################################################################################
    # types

    @property
    def types(self):
        # spath types
        return ([self.path_in_code] * len(self._path_in),
                [self.path_object_code] * len(self._path_object),
                [self.path_out_code] * len(self._path_out))

    @property
    def types_cont(self):
        typesc = []
        for t in self.types:
            typesc += t
        return typesc

    @property
    def gtypes(self):
        # generic types
        return self.types_in, self.types_object, self.types_out

    @property
    def gtypes_cont(self):
        gtypesc = []
        for t in self.gtypes:
            gtypesc += t
        return gtypesc

    @property
    @tupleify
    def etypes(self):
        # extended types
        for t, g in zip(self.types, self.gtypes):
            yield [t for t in zip(t, g)]

    @property
    def etypes_cont(self):
        return [t for t in zip(self.types_cont, self.gtypes_cont)]

    ####################################################################################################################

    @property
    def size(self):
        return sum(self.sizes)

    @property
    def sizes(self):
        return list(map(len, self._paths))

    @property
    def begins(self):
        return self.paths_cont[0]

    @property
    def ends(self):
        return self.paths_cont[-1]

    @property
    def has_in(self):
        return len(self._path_in) > 0

    @property
    def has_object(self):
        return len(self._path_object) > 0

    @property
    def has_out(self):
        return len(self._path_out) > 0

    ####################################################################################################################

    @tupleify
    def get_coords(self, smooth=None):
        # TODO: it is not used to get smooth coords but to get coords in general, conditionally smoothed
        # if smooth is not none applies smoothing
        if smooth is not None:
            for nr, coords in enumerate(self._make_smooth_coords(smooth)):
                if coords is None:
                    yield self.coords[nr]
                else:
                    yield coords
        else:
            for coords in self.coords:
                yield coords

    def _make_smooth_coords(self, smooth):
        return self.single_res_selection.coords_smooth(self._paths, smooth)

    def get_coords_cont(self, smooth=None):
        # returns coords as one array
        return make_default_array(np.vstack([c for c in self.get_coords(smooth) if len(c) > 0]))

    '''
    def apply_smoothing(self, smooth):
        # permament change!
        self.coords_in, self.coords_object, self.coords_out = self._make_smooth_coords(self.coords_cont, smooth)
    '''

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


class SinglePath(MacroMolPath):

    # __slots__ = "id _path_in _path_object _path_out _types_in _types_object _types_out _object_len single_res_selection".split()

    def is_single(self):
        return True

    def is_passing(self):
        return False


# TODO: there are problems with passing paths
class PassingPath(MacroMolPath):
    __slots__ = "_has_in_flag _has_out_flag".split()  # id _path_in _path_object _path_out _types_in _types_object _types_out _object_len single_res_selection".split()

    def __init__(self, path_id, paths, types):

        super(PassingPath, self).__init__(path_id, paths, types)

        self._has_in_flag = None
        self._has_out_flag = None

    def __getstate__(self):
        return self.id, self._path_in, self._path_object, self._path_out, self._types_in, self._types_object, self._types_out, self._object_len, self._has_in_flag, self._has_out_flag

    def __setstate__(self, state):
        self.id, self._path_in, self._path_object, self._path_out, self._types_in, self._types_object, self._types_out, self._object_len, self._has_in_flag, self._has_out_flag = state
        self.single_res_selection = SingleResidueSelection(self.id.id)

    def is_single(self):
        return False

    def is_passing(self):
        return True

    @property
    def has_in(self):
        return self._has_in_flag

    @has_in.setter
    def has_in(self, flag):
        self._has_in_flag = flag

    @property
    def has_out(self):
        return self._has_out_flag

    @has_out.setter
    def has_out(self, flag):
        self._has_out_flag = flag

    @property
    def coords_first_in(self):
        return self.coords_in[0]

    @property
    def paths_first_in(self):
        return self.path_in[0]

    @property
    def coords_last_out(self):
        return self.coords_in[-1]

    @property
    def paths_last_out(self):
        return self.path_in[-1]

    @property
    def types(self):
        # spath types
        return ([self.path_walk_code] * len(self._path_in),
                [],
                [])

    '''
    def get_coords(self, smooth=None):
        # TODO: it is not used to get smooth coords but to get coords in general, conditionally smoothed
        # if smooth is not none applies smoothing
        if smooth is not None:
            for nr,coords in enumerate(self._make_smooth_coords(smooth)):
                if coords is None:
                    yield self.coords[nr]
                else:
                    yield coords
        else:
            for coords in self.coords:
                yield coords

    '''

    def get_inlets(self):
        # Lets assume that if max/min_possible_frame is there is no inlet
        # how to check it?
        if self.has_in:
            yield Inlet(coords=self.coords[0][0],
                        inlet_type=(InletTypeCodes.surface, InletTypeCodes.incoming),
                        reference=self.id,
                        frame=self.paths_first_in)
        if self.has_out:
            yield Inlet(coords=self.coords[0][-1],
                        inlet_type=(InletTypeCodes.surface, InletTypeCodes.outgoing),
                        reference=self.id,
                        frame=self.paths_last_out)

    def remove_inlet(self, inlet_type):
        # only surface type can be removed
        if InletTypeCodes.surface == inlet_type[0]:
            if self.has_in and InletTypeCodes.incoming == inlet_type[1]:
                self.has_in = False
            elif self.has_out and InletTypeCodes.outgoing == inlet_type[1]:
                self.has_out = False
        # TODO: raise exceptions or at least warnings if wrong type of inlet is used


####################################################################################################################


class MasterPath(MacroMolPath):
    __slots__ = 'width_cont'.split()

    def __init__(self, sp, single_res_selection=None):
        super(MasterPath, self).__init__(sp.id, sp.paths, sp.gtypes)
        self.single_res_selection = single_res_selection
        self.width_cont = None

    def add_width(self, width):
        assert len(width) == self.size
        self.width_cont = width

    def __getstate__(self):
        return super(MasterPath, self).__getstate__() + (self.width_cont, self.single_res_selection)

    def __setstate__(self, state):
        self.id, self._path_in, self._path_object, self._path_out, self._types_in, self._types_object, self._types_out, self._object_len, self.width_cont, self.single_res_selection = state
