# import it as pmc?

import cPickle as pickle
import numpy as np
import os
import pymol
import tarfile
from pymol import cgo
from pymol import cmd

from aqueduct.geom import traces
from aqueduct.traj.paths import PathTypesCodes
from aqueduct.utils.helpers import list_blocks_to_slices
from aqueduct.visual.quickplot import cc, color_codes


# TODO: remove pymol form required dependencies - move it to pymol_real
# TODO: remove matplotlib from required dependencies

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


from aqueduct.utils.helpers import create_tmpfile


class SimpleTarWriteHelper(object):
    def __init__(self):
        self.tar_fh = None
        self.tmp_file = create_tmpfile()

    def open(self, filename):
        self.tar_fh = tarfile.open(filename, 'w:gz')

    def save_object2tar(self, obj, name):
        with open(self.tmp_file, 'w') as f:
            pickle.dump(obj, f)
        self.save_file2tar(self.tmp_file, name)

    def save_file2tar(self, filename, name):
        self.tar_fh.add(filename, arcname=name)

    def __del__(self):
        if self.tar_fh is not None:
            self.tar_fh.close()
        os.unlink(self.tmp_file)


class ConnectToPymol(object):
    cgo_line_width = 2.
    ct_pymol = 'pymol'
    ct_file = 'file'

    def __init__(self):
        self.connection_type = None  # possible types are pymol and file

        self.script_fh = None
        self.data_fh = SimpleTarWriteHelper()

    def init_pymol(self):
        pymol.finish_launching()
        cmd.set('cgo_line_width', ConnectToPymol.cgo_line_width)
        self.connection_type = self.ct_pymol

    def init_script(self, filename):
        self.script_fh = open(filename, 'w')
        data_filename = os.path.splitext(os.path.basename(filename))[0] + '.tar.gz'
        self.data_fh.open(data_filename)
        self.connection_type = self.ct_file

        # init lines, imports etc.
        self.script_fh.write('''print "Loading Aqua-Duct visualization..."
from pymol import cgo, cmd
cmd.set('cgo_line_width', %d)
from os import close, unlink
from os.path import splitext
import tarfile
import cPickle as pickle
from tempfile import mkstemp
fd, pdb_filename = mkstemp(suffix='pdb')
close(fd)
data_fh = tarfile.open("%s","r:gz")
def load_object(filename,name,state):
    print "Loading %s" % splitext(filename)[0]
    cmd.load_cgo(pickle.load(data_fh.extractfile(filename)), name, state)
    if state<2:
        cmd.refresh()
def load_pdb(filename,name,state):
    with open(pdb_filename,'w') as fpdb:
        fpdb.write(data_fh.extractfile(filename).read())
    cmd.load(pdb_filename, state=state, object=name)
''' % (self.cgo_line_width, data_filename,"%s","% s"))

    def add_cgo_object(self, name, cgo_object, state=None):
        if state is None:
            state = 1
        if self.connection_type == self.ct_pymol:
            cmd.load_cgo(cgo_object, str(name), state)
        elif self.connection_type == self.ct_file:
            obj_name = '%s_%d.dump' % (name, state)
            self.data_fh.save_object2tar(cgo_object, obj_name)

            self.script_fh.write('''load_object('%s', "%s", %d)''' % (obj_name, str(name), state))
            self.script_fh.write(os.linesep)
            # self.script_fh.write('''cmd.refresh()''')
            # self.script_fh.write(os.linesep)

    def del_cgo_object(self, name, state=None):
        raise NotImplementedError("This feature is not implemented yet.")

    def load_pdb(self, name, filename, state=None):
        if state is None:
            state = 1
        if self.connection_type == self.ct_pymol:
            cmd.load(filename, state=state, object=name)
        elif self.connection_type == self.ct_file:
            # save pdblile as string
            filename_new = '%s_%d.pdb' % (name, state)
            self.data_fh.save_file2tar(filename, filename_new)
            self.script_fh.write('''load_pdb("%s","%s",%d)''' % (filename_new,name,state))
            self.script_fh.write(os.linesep)

    def orient_on(self, name):
        if self.connection_type == self.ct_pymol:
            cmd.orient(name)
        elif self.connection_type == self.ct_file:
            self.script_fh.write('''cmd.orient("%s")''' % name)
            self.script_fh.write(os.linesep)

    def __del__(self):
        if self.connection_type == self.ct_file:
            self.script_fh.write('''unlink(pdb_filename)
print "Aqua-Duct visualization loaded."''')
            self.script_fh.write(os.linesep)
            self.script_fh.close()


class SinglePathPlotter(object):
    def __init__(self, pymol_connector, linearize=None):

        self.cgo_lines = BasicPymolCGOLines()
        self.cgo_spheres = BasicPymolCGOSpheres()
        self.cgo_pointers = BasicPymolCGOPointers()

        self.linearize = linearize

        self.pymol_connector = pymol_connector

    # TODO: take care of proper colors handling for smoothed and not smoothed traces!

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

        # create slices
        if smooth:
            sls = tuple(list_blocks_to_slices(spath.types_cont))
        else:
            sls = tuple(list_blocks_to_slices(spath.etypes_cont))

        # create traces
        traces_list = tuple([coords_cont[sl] for sl in sls])
        new_line = False
        for trace, sl in zip(traces.midpoints(traces_list), sls):
            if len(trace) > 0:
                # now trace has midpoints
                # get type and etype
                t = spath.types_cont[sl][0]
                et = spath.etypes_cont[sl][0]
                if smooth:
                    et = et[0]
                # plot, if allowed
                if (plot_in and t == PathTypesCodes.path_in_code) or (
                            plot_object and t == PathTypesCodes.path_object_code) or (
                            plot_out and t == PathTypesCodes.path_out_code):
                    # get color
                    c = color_codes(et)
                    # now, it is possible to linearize!
                    if smooth and self.linearize:
                        trace = traces.LinearizeRecursiveVector(self.linearize)(trace)
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
        self.pymol_connector.add_cgo_object(name, self.cgo_lines.get(), state=state)

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
        self.pymol_connector.add_cgo_object(name, self.cgo_pointers.get(), state=state)

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

        self.pymol_connector.add_cgo_object(name, self.cgo_spheres.get(), state=state)

    def convexhull(self, chull,
                   color='m',
                   name='convexhull',
                   state=None):

        if isinstance(color, str):
            color = cc(color)
        if state is None:
            state = 1

        self.cgo_lines.clean()
        for nr, facet in enumerate(chull.facets):
            if nr > 0:
                self.cgo_lines.new()
            self.cgo_lines.add(facet, color=color)

        self.pymol_connector.add_cgo_object(name, self.cgo_lines.get(), state=state)
