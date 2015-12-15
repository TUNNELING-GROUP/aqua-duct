'''
Created on Dec 10, 2015

@author: tljm
'''

from itertools import izip

def is_iterable(l):
    try:
        _ = (e for e in l)
        return True
    except TypeError:
        pass
    return False


def listify(gen):
    # http://argandgahandapandpa.wordpress.com/2009/03/29/python-generator-to-list-decorator/
    # improved by tljm: for non iterable objects it returns list of one element: the object
    def patched(*args, **kwargs):
        obj = gen(*args, **kwargs)
        if is_iterable(obj):
            return list(obj)
        return [obj]

    return patched

def sortify(gen):
    # http://argandgahandapandpa.wordpress.com/2009/03/29/python-generator-to-list-decorator/
    # improved by tljm: for non iterable objects it returns list of one element: the object
    def patched(*args, **kwargs):
        obj = gen(*args, **kwargs)
        if is_iterable(obj):
            obj = list(obj)
            obj.sort()
            return obj
        return [obj]

    return patched

   
# paths/list manipulations

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

#TODO: create utils module

class GenericPath(object):
    # object to store paths... is it required?
    # could be nice to have as it can be armed with methods for getting in and out paths
    object_name = 'o'
    scope_name = 's'
    out_name = 'n'

    def __init__(self):
        
        self.types = []
        self.frames = []
        self.coords = []
        
    def add_coord(self,coord):
        self.coords.append(coord)
    
    def add_object(self,frame):
        self.types.append(self.object_name)
        self.frames.append(frame)
        
    def add_scope(self,frame):
        self.types.append(self.scope_name)
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


    def find_paths(self):
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
            yield path_in,path_core,path_out
        path_in = []
        path_core = []
        for path_out in paths_out:
            yield path_in,path_core,path_out

    def get_single_path_coords(self,spath,smooth=None):
        # returns coordinates for single path
        # single path coprises of in,scope,out parts
        # smooth has to be callable accepting coords
        
        p_in,p_scope,p_out = spath
        
        # get smoothed coords
        if smooth:
            coords = smooth(self.coords).tolist()
        else:
            coords = self.coords
            
        if len(self.frames) == len(coords):
            # not full trajectory
            in_ = []
            for f in p_in:
                in_.append(coords[self.frames.index(f)])
            scope_ = []
            for f in p_scope:
                scope_.append(coords[self.frames.index(f)])
            out_ = []
            for f in p_out:
                out_.append(coords[self.frames.index(f)])
        else:
            # full trajectory
            in_ = []
            for f in p_in:
                in_.append(coords[f])
            scope_ = []
            for f in p_scope:
                scope_.append(coords[f])
            out_ = []
            for f in p_out:
                out_.append(coords[f])

        return in_,scope_,out_
            
if __name__ == "__main__":
    
    p = GenericPath()
    p.types = ['s', 'o', 's', 'o', 'o', 'o', 'o', 's', 's', 's', 's', 's', 's']
    p.frames = [0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
    
    for t,f in izip(p.types,p.frames):
        print t,f
    
    print list(p.get_paths_in())
    print list(p.get_paths_out())
    print list(p.find_paths())