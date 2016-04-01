# import it as pmc?

import numpy as np

from aqueduct.geom import traces
from aqueduct.visual.quickplot import cc
from aqueduct.traj.paths import PathTypesCodes, GenericPathTypeCodes
from aqueduct.utils.helpers import list_blocks_to_slices

import pymol

from pymol import cgo
from pymol import cmd


class BasicPymolCGO(object):
    cgo_entity_begin = []
    cgo_entity_end = []

    def __init__(self):
        self.cgo_entity = None
        self.previous = None
        self.clean()

    def clean(self):
        self.previous = None
        self.cgo_entity = []
        self.cgo_entity.extend(self.cgo_entity_begin)

    def new(self):
        self.previous = None

    def get(self):
        return self.cgo_entity + self.cgo_entity_end


class BasicPymolCGOLines(BasicPymolCGO):
    cgo_entity_begin = [cgo.BEGIN, cgo.LINES]
    cgo_entity_end = [cgo.END]

    def add(self, coords=None, color=None):

        if color is not None:
            self.cgo_entity.append(cgo.COLOR)
            self.cgo_entity.extend(map(float, color))

        if coords is not None:
            for nr, coord in enumerate(coords):
                if self.previous is not None:
                    self.cgo_entity.append(cgo.VERTEX)
                    self.cgo_entity.extend(map(float, self.previous))

                    self.cgo_entity.append(cgo.VERTEX)
                    self.cgo_entity.extend(map(float, coord))
                self.previous = coord


class BasicPymolCGOSpheres(BasicPymolCGO):
    cgo_entity_begin = []
    cgo_entity_end = []

    def add(self, coords=None, radius=None, color=None):
        # color to colors...
        if color is not None:
            color = np.matrix(color).A
        if radius is not None:
            radius = np.matrix(radius).A1

        if coords is not None:
            for nr, coord in enumerate(coords):
                if color is not None:
                    if len(color) > 1:
                        c = color[nr]
                    else:
                        c = color[0]
                    self.cgo_entity.append(cgo.COLOR)
                    self.cgo_entity.extend(map(float, c))
                self.cgo_entity.append(cgo.SPHERE)
                self.cgo_entity.extend(map(float, coord))
                if radius is not None:
                    if len(radius) > 1:
                        r = radius[nr]
                    else:
                        r = radius
                else:
                    r = 1.
                self.cgo_entity.append(float(r))


class ConnectToPymol(object):
    cgo_line_width = 2.

    @staticmethod
    def init_pymol():
        pymol.finish_launching()
        cmd.set('cgo_line_width', ConnectToPymol.cgo_line_width)

    @staticmethod
    def add_cgo_object(name, cgo_object, state=None):
        if state is None:
            state = 1
        cmd.load_cgo(cgo_object, str(name), state)

    @staticmethod
    def del_cgo_object(name, state=None):
        raise NotImplementedError("This feature is not implemented yet.")

    @staticmethod
    def load_pdb(name, filename, state=None):
        if state is None:
            state = 1
        cmd.load(filename, state=state, object=name)


class SinglePathPlotter(object, PathTypesCodes):
    def __init__(self):

        self.cgo_lines = BasicPymolCGOLines()
        self.cgo_spheres = BasicPymolCGOSpheres()

    def add_single_path_continous_trace(self,
                                        spath,
                                        smooth=None,
                                        color=('r', 'g', 'b', 'y'),
                                        plot_in=True,
                                        plot_object=True,
                                        plot_out=True,
                                        **kwargs):

        self.cgo_lines.new()

        for nr, trace in enumerate(traces.midpoints(spath.get_smooth_coords(smooth))):
            # mid points!
            if len(trace) > 0:
                if (nr == 0 and plot_in) or (nr == 1 and plot_object) or (nr == 2 and plot_out):
                    if nr == 1:
                        # this is special case of object
                        sts = list(list_blocks_to_slices(spath.types_object))
                        for strace, gtype in zip(traces.midpoints(tuple([trace[s, :] for s in sts])),
                                                 [spath.types_object[s] for s in sts]):
                            if gtype[0] == self.path_object_code:
                                c = color[1]
                            else:
                                c = color[3]
                            self.cgo_lines.add(strace, cc(c))
                    else:
                        self.cgo_lines.add(trace, cc(color[nr]))
                else:
                    self.cgo_lines.new()
            else:
                self.cgo_lines.new()

    def paths_trace(self,
                    spaths,
                    smooth=None,
                    color=('r', 'g', 'b', 'y'),
                    name='paths',
                    state_function=None,
                    **kwargs):

        # if state_function is None
        if state_function is None:
            self.cgo_lines.clean()
            for nr, spath in enumerate(spaths):
                self.add_single_path_continous_trace(spath, smooth=smooth, color=color, **kwargs)
            ConnectToPymol.add_cgo_object(name, self.cgo_lines.get(), state=1)
        # else if state_function is not None
        else:
            # first lets find min and max state
            min_state = []
            max_state = []
            for nr, spath in enumerate(spaths):
                mins, maxs = state_function(nr, spath)
                min_state.append(mins)
                max_state.append(maxs)
            min_state = min(min_state)
            max_state = max(max_state)
            # now, loop over possible states and find paths that fits into it
            for state in range(min_state, max_state + 1):
                self.cgo_lines.clean()
                # now find paths that fits
                for nr, spath in enumerate(spaths):
                    mins, maxs = state_function(nr, spath)
                    if state >= mins and state <= maxs:
                        self.add_single_path_continous_trace(spath, smooth=smooth, color=color, **kwargs)
                ConnectToPymol.add_cgo_object(name, self.cgo_lines.get(), state=state)

    def scatter(self, coords, radius=0.4, color='r', name='scatter', state=None):

        if isinstance(color, str):
            color = cc(color)

        if state is None:
            state = 1

        self.cgo_spheres.clean()

        self.cgo_spheres.add(coords=coords, radius=radius, color=color)

        ConnectToPymol.add_cgo_object(name, self.cgo_spheres.get(), state=state)
