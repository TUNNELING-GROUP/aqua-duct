'''
Created on Dec 10, 2015

@author: tljm
'''

from itertools import izip
import numpy as np
from aqueduct.utils.helpers import tupleify, sortify, is_iterable, listify

   
########################################################################################################################
# paths/list manipulations
# following part of code comes directly from the very initial tcl/vmd implementation
# all following functions should return list

def union(a,b):
    return [aa for aa in a if aa in b]

def glue(a,b):
    if a[-1] >= b[0]:
        g = list(set(a+b))
        g.sort()
        return g
    return []

@listify
def xor(a,b):
    ab = union(a,b)
    for e in glue(a,b):
        if e not in ab:
            yield e

def left(a,b):
    return union(a,xor(a,b))

def right(a,b):
    return union(b,xor(a,b))

########################################################################################################################

class PathTypesCodes():
    path_in_code = 'i'
    path_object_code = 'c'
    path_out_code = 'o'
    

class GenericPaths(object):
    # object to store paths... is it required?
    # could be nice to have as it can be armed with methods for getting in and out paths
    object_name = 'c'
    scope_name = 's'
    out_name = 'n'

    def __init__(self,id,min_pf=None,max_pf=None):

        # id is any type of object; it is used as identifier

        self.id = id

        self.types = []
        self.frames = []
        self.coords = []

        # following is required to correct in and out paths that begin or end in scope and
        # begin or end at the very begining of MD or at very end of MD

        self.max_possible_frame = max_pf
        self.min_possible_frame = min_pf

    def add_coord(self,coord):
        self.coords.append(coord)
    
    def add_object(self,frame):
        self.add_type(frame,self.object_name)

    def add_scope(self,frame):
        self.add_type(frame,self.scope_name)

    def add_type(self,frame,ftype):
        self.types.append(ftype)
        self.frames.append(frame)

    @property
    def max_frame(self):
        return max(self.frames)

    @property
    def min_frame(self):
        return min(self.frames)
        
    #xrange(self.max_frame,self.min_frame-1,-1)
        
    def get_paths_in(self):
        frames_range = xrange(self.max_frame,self.min_frame-1,-1)
        for path in self.get_paths_for_frames_range(frames_range):
            path.sort()
            yield path
    
    def get_paths_out(self):
        frames_range = xrange(self.min_frame,self.max_frame+1,1)
        for path in self.get_paths_for_frames_range(frames_range):
            path.sort()
            yield path
    
    @sortify    
    def get_paths_for_frames_range(self,frames_range):
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


    def find_paths(self,fullonly=False):
        paths_out = list(self.get_paths_out())
        paths_in = list(self.get_paths_in())
        
        for path_in in paths_in:
            path_out = paths_out[0]
            path_glue = glue(path_in,path_out)
            if len(path_glue) > 0:
                path_out = paths_out.pop(0)
                path_core = union(path_in,path_out)
                path_in = left(path_in,path_out)
                path_out = right(path_core,path_out)
            else:
                path_core = []
                path_out = []
            if len(path_in) == 0 or len(path_core) == 0 or len(path_out) == 0:
                if fullonly:
                    continue
            #TODO: make separate method/function for that
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
            yield path_in,path_core,path_out
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
                yield path_in,path_core,path_out
            
    def find_paths_coords(self,fullonly=False):
        for path in self.find_paths(fullonly=fullonly):
            yield self.get_single_path_coords(path)

    def find_paths_types(self,fullonly=False):
        for path in self.find_paths(fullonly=fullonly):
            yield self.get_single_path_types(path)

    def get_single_path_coords(self,spath):
        # returns coordinates for single path
        # single path comprises of in,scope,out parts

        p_in,p_object,p_out = spath
        
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

        return in_,object_,out_

    def get_single_path_types(self,spath):
        # returns typess for single path
        # single path comprises of in,scope,out parts

        #TODO: join it with get_single_path_coords

        p_in,p_object,p_out = spath

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

        return in_,object_,out_


def yield_single_paths(gps, fullonly=False, progress=False):
    # iterates over gps - list of GenericPaths objects and transforms them in to SinglePath objects
    for nr,gp in enumerate(gps):
        id = gp.id
        for paths,coords,types in zip(gp.find_paths(fullonly=fullonly),
                                      gp.find_paths_coords(fullonly=fullonly),
                                      gp.find_paths_types(fullonly=fullonly)):
            if progress:
                yield SinglePath(id,paths,coords,types),nr

            else:
                yield SinglePath(id,paths,coords,types)

def list_blocks_to_slices(l):
    n = len(l)
    if n in [0,1]:
        yield slice(None,None,None)
    if n > 1:
        prev = l[0]
        prev_nr = 0
        for nr,e in enumerate(l[1:]):
            if e == prev:
                continue
            yield slice(prev_nr,nr+1,1)
            prev = e
            prev_nr = nr+1
        yield slice(prev_nr,nr+2,1)


class InletTypeCodes:
    inlet_in_code = 'inin'
    inlet_out_code = 'inout'


class SinglePath(object,PathTypesCodes,InletTypeCodes):
    # special class
    # represents one path


    empty_coords = np.zeros((0,3))

    def __init__(self, id, paths, coords, types):

        self.id = id
        self.path_in,self.path_object,self.path_out = paths
        self.coords_in,self.coords_object,self.coords_out = map(np.array,coords)
        self.types_in, self.types_object, self.types_out = types

        #return np.vstack([c for c in self._coords if len(c) > 0])

    @property
    def size(self):
        return sum(map(len,self.paths))

    @property
    def coords(self):
        return self.coords_in,self.coords_object,self.coords_out

    @property
    def coords_cont(self):
        # returns coords as one array
        return np.vstack([c for c in self.coords if len(c) > 0])


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
        return [(inlet,{0:self.inlet_in_code,1:self.inlet_out_code}[nr]) for nr,inlet in enumerate((self.coords_first_in,self.coords_last_out)) if inlet is not None]

    @property
    def paths(self):
        return self.path_in,self.path_object,self.path_out

    @property
    def paths_cont(self):
        return self.path_in+self.path_object+self.path_out

    @property
    def types_cont(self):
        return ([self.path_in_code]*len(self.path_in))+([self.path_object_code]*len(self.path_object))+([self.path_out_code]*len(self.path_out))

    @property
    def gtypes(self):
        return self.types_in, self.types_object, self.types_out

    @property
    def gtypes_cont(self):
        return self.types_in + self.types_object + self.types_out

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

    @tupleify
    def get_smooth_coords(self,smooth):
        # smooth should be callable and should return an object of length equal to submitted one
        # get continuous coords
        if smooth:
            coords_smooth = smooth(self.coords_cont)
        else:
            coords_smooth = self.coords_cont
        # now lets return tupple of coords
        nr = 0
        for path in self.paths:
            if len(path) > 0:
                yield np.array(coords_smooth[nr:nr+len(path)])
                nr += len(path)
            else:
                yield self.empty_coords
    
    def apply_smoothing(self,smooth):
        # applies smoothing
        # this is permamant change
        self.coords_in,self.coords_object,self.coords_out = self.get_smooth_coords(smooth)
        


if __name__ == "__main__":

    import numpy as np


    print "in starts at 0"
    gp = GenericPaths(id=0,min_pf=0,max_pf=100)

    frames = range(20)
    types = gp.scope_name * 10 + gp.object_name * 10
    coords = np.random.randn(20,3)
    for f,t,c in zip(frames,types,coords):
        gp.add_coord(c)
        gp.add_type(f,t)

    print list(gp.find_paths())
    print list(gp.find_paths_types())


    print "in starts at 1"
    gp = GenericPaths(id=0,min_pf=0,max_pf=100)

    frames = range(1,21)
    types = gp.scope_name * 10 + gp.object_name * 10
    coords = np.random.randn(20,3)
    for f,t,c in zip(frames,types,coords):
        gp.add_coord(c)
        gp.add_type(f,t)

    print list(gp.find_paths())
    print list(gp.find_paths_types())

    print "out ends at 29"
    gp = GenericPaths(id=0,min_pf=0,max_pf=10)

    frames = range(0,30)
    types = gp.scope_name * 10 + gp.object_name * 10 + gp.scope_name * 10
    coords = np.random.randn(30,3)
    for f,t,c in zip(frames,types,coords):
        gp.add_coord(c)
        gp.add_type(f,t)

    print list(gp.find_paths())
    print list(gp.find_paths_types())
