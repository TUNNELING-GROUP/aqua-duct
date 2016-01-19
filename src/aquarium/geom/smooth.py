'''
Created on Dec 15, 2015

@author: tljm
'''

from aquarium.utils.helpers import arrayify

import numpy as np

class Smooth(object):

    def __init__(self,recursive=None):

        self.recursive = recursive

    def smooth(self,coords):
        raise NotImplementedError("Missing implementation")

    def __call__(self,coords):
        if self.recursive:
            coords_smooth = None
            for r in xrange(self.recursive):
                if coords_smooth is not None:
                    coords_smooth = self.smooth(coords_smooth)
                else:
                    coords_smooth = self.smooth(coords)
            return coords_smooth
        return self.smooth(coords)

class WindowSmooth(Smooth):
    
    def __init__(self,window=5,function=np.mean,**kwargs):
        Smooth.__init__(self,**kwargs)
        self.window = window
        self.function = function

    @arrayify
    def smooth(self,coords):
        
        n = len(coords)
        
        for pos in xrange(n):
            if pos < self.window:
                lo = None
            else:
                lo = pos-self.window
            if n - 1 - pos < self.window:
                hi = None
            else:
                hi = pos+self.window
            
            yield self.function(coords[slice(lo,hi,None)],0).tolist()
            