
# import it as pmc?

from aquarium.geom import traces
from aquarium.visual.quickplot import cc

import pymol

from pymol import cgo
from pymol import cmd

class BasicPymolCGOLines(object):

    def __init__(self):

        self.line = None
        self.previous = None
        self.clean()

    def clean(self):
        self.previous = None
        self.line = [cgo.BEGIN,cgo.LINES]

    def new(self):
        self.previous = None


    def add(self,coords=None,color=None):

        if color is not None:
            self.line.append(cgo.COLOR)
            self.line.extend(map(float,color))

        if coords is not None:
            for nr,coord in enumerate(coords):
                if self.previous is not None:
                    self.line.append(cgo.VERTEX)
                    self.line.extend(map(float,self.previous))

                    self.line.append(cgo.VERTEX)
                    self.line.extend(map(float,coord))
                self.previous = coord


    def get(self):
        return self.line+[cgo.END]



class ConnectToPymol(object):

    @staticmethod
    def init_pymol():
        pymol.finish_launching()

    @staticmethod
    def add_cgo_object(name,cgo_object,state=None):
        if state is None:
            state = 1
        cmd.load_cgo(cgo_object,str(name), state)

    @staticmethod
    def del_cgo_object(name,state=None):

        raise NotImplementedError("This feature is not implemented yet.")



class SinglePathPlotter(object):

    def __init__(self):

        self.cgo_object = BasicPymolCGOLines()


    def add_single_path_continous_trace(self,
                                       spath,
                                       smooth=None,
                                       color=('r','g','b'),
                                       plot_in=True,
                                       plot_object=True,
                                       plot_out=True,
                                       **kwargs):

        self.cgo_object.new()

        for nr,trace in enumerate(traces.midpoints(spath.get_smooth_coords(smooth))):
            # mid points!
            if len(trace) > 0:
                if (nr == 0 and plot_in) or (nr == 1 and plot_object) or (nr == 2 and plot_out):
                    self.cgo_object.add(trace,cc(color[nr]))
                else:
                    self.cgo_object.new()
            else:
                self.cgo_object.new()



    def paths_trace(self,
                    spaths,
                    smooth=None,
                    color=('r','g','b'),
                    name='paths',
                    state_function=None,
                    **kwargs):

        # if state_function is None
        if state_function is None:
            self.cgo_object.clean()
            for nr,spath in enumerate(spaths):
                self.add_single_path_continous_trace(spath,smooth=smooth,color=color,**kwargs)
            ConnectToPymol.add_cgo_object(name,self.cgo_object.get(),state=1)
        # else if state_function is not None
        else:
            # first lets find min and max state
            min_state = []
            max_state = []
            for nr,spath in enumerate(spaths):
                mins,maxs = state_function(nr,spath)
                min_state.append(mins)
                max_state.append(maxs)
            min_state = min(min_state)
            max_state = max(max_state)
            # now, loop over possible states and find paths that fits into it
            for state in range(min_state,max_state+1):
                self.cgo_object.clean()
                # now find paths that fits
                for nr,spath in enumerate(spaths):
                    mins,maxs = state_function(nr,spath)
                    if state >= mins and state <= maxs:
                        self.add_single_path_continous_trace(spath,smooth=smooth,color=color,**kwargs)
                ConnectToPymol.add_cgo_object(name,self.cgo_object.get(),state=state)



