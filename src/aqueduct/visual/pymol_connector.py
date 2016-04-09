# import it as pmc?

import pymol
from pymol import cgo
from pymol import cmd

import numpy as np

from aqueduct.geom import traces
from aqueduct.traj.paths import PathTypesCodes
from aqueduct.utils.helpers import list_blocks_to_slices
from aqueduct.visual.quickplot import cc, color_codes


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


# [cgo.CONE] + xyz3 + xyz2 + [hradius, 0.0] + color2 + color2

class BasicPymolCGOPointers(BasicPymolCGO):
    cgo_entity_begin = []
    cgo_entity_end = []

    def add_cone(self, coords1=None, coords2=None, radius1=None, radius2=None, color1=None, color2=None):
        # color to colors... ???
        if coords1 is not None and coords2 is not None:
            self.cgo_entity.append(cgo.CONE)
            self.cgo_entity.extend(map(float, coords1))
            self.cgo_entity.extend(map(float, coords2))
            self.cgo_entity.append(float(radius1))
            self.cgo_entity.append(float(radius2))
            self.cgo_entity.extend(map(float, color1))
            self.cgo_entity.extend(map(float, color2))

            self.cgo_entity.extend([cgo.NULL, cgo.POINTS])

    def add_pointer(self, point=None, direction=None, length=None, color=None, reverse=False):
        vec = point - direction
        vec_len = np.sqrt(np.sum((vec) ** 2))
        vec = (vec / vec_len) * length
        vec = point + vec
        if reverse:
            self.add_cone(coords1=point, coords2=vec, radius1=length / 3., radius2=0, color1=color, color2=color)
        else:
            self.add_cone(coords1=vec, coords2=point, radius1=length / 3., radius2=0, color1=color, color2=color)


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


class SinglePathPlotter(object):
    def __init__(self):

        self.cgo_lines = BasicPymolCGOLines()
        self.cgo_spheres = BasicPymolCGOSpheres()
        self.cgo_pointers = BasicPymolCGOPointers()

    def add_single_path_continous_trace(self,
                                        spath,
                                        smooth=None,
                                        plot_in=True,
                                        plot_object=True,
                                        plot_out=True,
                                        **kwargs):

        self.cgo_lines.new()

        # get coords
        coords_cont = spath.get_coords_cont(smooth)
        # create traces
        traces_list = tuple([coords_cont[sl] for sl in list_blocks_to_slices(spath.etypes_cont)])

        new_line = False
        for trace, sl in zip(traces.midpoints(traces_list), list_blocks_to_slices(spath.etypes_cont)):
            if len(trace) > 0:
                # now trace has midpoints
                # get type and etype
                t = spath.types_cont[sl][0]
                et = spath.etypes_cont[sl][0]
                # plot, if allowed
                if (plot_in and t == PathTypesCodes.path_in_code) or (
                    plot_object and t == PathTypesCodes.path_object_code) or (
                    plot_out and t == PathTypesCodes.path_out_code):
                    # get color
                    c = color_codes(et)
                    self.cgo_lines.add(trace, cc(c))
                    # new_line = True

            if new_line:
                self.cgo_lines.new()
                new_line = False

    def paths_trace(self, spaths,
                    smooth=None,
                    name='paths',
                    state=None,
                    **kwargs):

        # if state_function is None
        if state is None:
            state = 1
        self.cgo_lines.clean()
        for nr, spath in enumerate(spaths):
            self.add_single_path_continous_trace(spath, smooth=smooth, **kwargs)
        ConnectToPymol.add_cgo_object(name, self.cgo_lines.get(), state=state)

    def paths_inlets(self, spaths,
                     smooth=None,
                     color=None,
                     plot_in=True,
                     plot_out=True,
                     name='in-out-let',
                     state=None,
                     **kwargs):

        if state is None:
            state = 1
        self.cgo_pointers.clean()
        for nr, spath in enumerate(spaths):
            coords = spath.get_coords_cont(smooth=smooth)
            etypes = spath.etypes_cont
            if plot_in:
                if color is None:
                    c = color_codes(etypes[0])
                else:
                    c = color
                self.cgo_pointers.add_pointer(point=coords[0],
                                              direction=coords[1],
                                              length=1.,
                                              color=cc(c))
            if plot_out:
                if color is None:
                    c = color_codes(etypes[-1])
                else:
                    c = color
                self.cgo_pointers.add_pointer(point=coords[-1],
                                              direction=coords[-2],
                                              length=1.,
                                              color=cc(c),
                                              reverse=True)
        ConnectToPymol.add_cgo_object(name, self.cgo_pointers.get(), state=state)

    def scatter(self, coords,
                radius=0.4,
                color='r',
                name='scatter',
                state=None):

        if isinstance(color, str):
            color = cc(color)

        if state is None:
            state = 1

        self.cgo_spheres.clean()

        self.cgo_spheres.add(coords=coords, radius=radius, color=color)

        ConnectToPymol.add_cgo_object(name, self.cgo_spheres.get(), state=state)
