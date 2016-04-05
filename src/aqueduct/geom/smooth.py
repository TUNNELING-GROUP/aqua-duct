'''
Created on Dec 15, 2015

@author: tljm
'''

from aqueduct.geom import traces
from aqueduct.utils.helpers import arrayify

import numpy as np

class Smooth(object):

    def __init__(self,recursive=None,**kwargs):

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
                hi = pos + 1
            else:
                lo = pos-self.window

                if n - 1 - pos < self.window:
                    lo = pos
                    hi = None
                else:
                    hi = pos+self.window

            yield self.function(coords[slice(lo,hi,None)],0).tolist()

class MaxStepSmooth(Smooth):

    def __init__(self, sigma=2, **kwargs):
        Smooth.__init__(self, **kwargs)
        self.sigma = sigma

    @arrayify
    def smooth(self, coords):

        n = len(coords)
        cdiff = traces.diff(coords)

        step = np.std(cdiff)*self.sigma

        current_step = 0
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
                        for coord in traces.tracepoints(last_coord,current_coord,to_yield_count):
                            yield coord
                    yield current_coord
                    last_coord = current_coord
                    to_yield_count = 0
                else:
                    current_step = traces.diff(np.vstack((current_coord,last_coord)))
                    if current_step > step:
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

    def __init__(self,**kwargs):
        Smooth.__init__(self, **kwargs)

        self.window = WindowSmooth(**kwargs)
        self.mss = MaxStepSmooth(**kwargs)

    def smooth(self,coords):
        return self.window(self.mss(coords))
