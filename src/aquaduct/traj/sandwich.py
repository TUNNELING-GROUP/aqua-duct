# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2018-2019  Tomasz Magdziarz, Michał Banas <info@aquaduct.pl>
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

from aquaduct import logger
from aquaduct.utils.helpers import version_parser
import MDAnalysis as mda

################################################################################
# Check MDAnalysis version

def mda_ver():
    return version_parser(mda.__version__)


# FIXME: do it according to user options
if mda_ver() < version_parser('0.16'):
    logger.error('Unsupported MDAnalysis version %s; should be 0.16.2 or > 0.19.', mda.__version__)
    raise NotImplementedError('Unsupported MDAnalysis version %s; should be 0.16.2 or > 0.19.' % mda.__version__)

if mda_ver() < version_parser('0.17') and mda_ver() < version_parser('0.20'):
    logger.warning('Unsupported MDAnalysis version %s.', mda.__version__)
    logger.warning('Trying to mitigate potential problems by setting `use_periodic_selections = False`.')
    mda.core.flags['use_periodic_selections'] = False

################################################################################
# rest of imports

import re
from os.path import splitext
from os import pathsep
from collections import OrderedDict, namedtuple

import numpy as np

from MDAnalysis.topology.core import guess_atom_element

from aquaduct.utils.helpers import SmartRange, SmartRangeIncrement
from aquaduct.geom.convexhull import SciPyConvexHull, are_points_within_convexhull
from aquaduct.utils.helpers import arrayify, create_tmpfile, tupleify
from aquaduct.utils.maths import make_default_array
from aquaduct.apps.data import GCS, CRIC
from aquaduct.utils.maths import defaults

################################################################################
# import or create memory decorator

if GCS.cachedir:
    from joblib import Memory
    memory_cache = Memory(cachedir=GCS.cachedir,
                          verbose=0)
    # mmap have to be switched off, otherwise smoothing does not work properly
    # memory_cache = Memory(cachedir=GCS.cachedir, mmap_mode='r', verbose=0)
    memory = memory_cache.cache
elif GCS.cachemem:
    from aquaduct.utils.helpers import memory_in_memory as memory
else:
    from aquaduct.utils.helpers import noaction as memory


################################################################################
# trajectory window object

class Window(object):
    def __init__(self, start, stop, step):
        self.start = self._none_or_int(start)
        self.stop = self._none_or_int(stop)
        self.step = self._none_or_int(step)

    @staticmethod
    def _none_or_int(nr):
        if nr is not None:
            return int(nr)

    def __repr__(self):
        return "Window(%r:%r:%r)" % (self.start, self.stop, self.step)

    def correct(self, real_frame_no):
        if self.start is None:
            self.start = 0
        elif self.start < 0:
            self.start = 0
        elif self.start > real_frame_no:
            self.start = real_frame_no - 1

        if self.stop is None:
            self.stop = real_frame_no - 1
        elif self.stop < 0:
            self.stop = 0
        elif self.stop > real_frame_no:
            self.stop = real_frame_no - 1

        if self.step is None:
            self.step = 1
        elif self.step < 0:
            self.step = 1

    def range(self, reverse=False):
        # returns range object
        if reverse:
            return range(self.stop, self.start - 1, -1 * self.step)
        return range(self.start, self.stop + 1, self.step)

    def get_real(self, frame):
        # recalculates real frame
        return self.start + frame * self.step

    def len(self):
        # lenght of window
        return len(self.range())

    def split(self, slices=None):
        assert slices > 0
        N = int(max(1, np.floor(self.len() / float(slices))))
        for start in list(self.range())[::(N + 1)]:
            yield Window(start, min(self.stop, start + (N) * self.step), self.step)


################################################################################
# MasterReader

# engine problem
# Two options are currently available:
# MDAnalysis
# MDTraj

class OpenReaderTraj(namedtuple('OpenReaderTraj', 'topology trajectory number window engine')):
    def open(self):
        # convenience function
        return open_traj_reader(self)


class MasterReader(object):

    # only one MasterReader object can be used
    # it does not use ANY direct call to ANY MD access software

    open_reader_traj = {}

    topology = ['']
    trajectory = ['']
    window = None  # this is full window

    sandwich_mode = None
    engine_name = 'mda'
    threads = 1
    threads_multiply = 1
    edges = []

    def __del__(self):
        self.reset()

    def __call__(self, topology, trajectory, window=None, sandwich=False, threads=1):
        """
        :param list topology:  List of topologies. Each element is a file name.
        :param list trajectory: List of trajectories. Each element is a file name.
        :param Window window: Frames window to read.
        :param bool sandwich: Flag for setting sandwich mode.

        If no sandiwch mode is used, number of topologies has to be precisely 1.
        In sandwich mode it can be either 1 or equal to the number of trajectory files.
        """

        if not isinstance(topology, list):
            topology = [t.strip() for t in topology.split(pathsep)]
        self.topology = topology
        if not isinstance(trajectory, list):
            trajectory = [t.strip() for t in trajectory.split(pathsep)]
        self.trajectory = trajectory

        self.window = window
        self.sandwich_mode = sandwich

        if len(self.topology) > 1:
            assert self.sandwich_mode, "Multiple topologies possible in sandwich mode only."
            assert len(self.topology) == len(
                self.trajectory), "Number of topologies must be 1 or be equal to number of trajectories."

        self.window.correct(self.real_number_of_frames())  # this corrects window
        self.reset()  # assert window correction clear all opened trajs

        self.threads = threads

    def reset(self):
        self.open_reader_traj = {}

    '''
    def __getstate__(self):
        # if pickle dump is required, this will not be used in the future
        # do not pass open_reader_traj
        return dict(((k, v) for k, v in self.__dict__.items() if k not in ['open_reader_traj']))

    def __setstate__(self, state):
        # if pickle dump is required, this will not be used in the future
        self.__dict__ = state
        self.open_reader_traj = {}

    '''
    def getrecallstate(self):
        # TODO: to be removed?
        return dict(topology=self.topology,
                    trajectory=self.trajectory,
                    window=self.window,
                    sandwich=self.sandwich_mode,
                    threads=self.threads)

    '''
    @property
    def engine(self):
        # returns engine used for accessing trajectory
        # this returns object that inherits from ReaderTraj
        if self.engine_name == 'mda':
            return ReaderTrajViaMDA
        raise NotImplementedError
    '''

    def engine(self, topology, trajectory, number, window):
        '''
        :param str topology: Topology file name.
        :param list trajectory: List of trajectories. Each element is a file name. Alternatively, trajectory file name.
        :param int number: Number of the layer.
        :param Window window: Frames window.
        :return:
        '''
        assert self.engine_name in ['mda']
        return OpenReaderTraj(topology, trajectory, number, window, self.engine_name)

    def __repr__(self):
        sandwich = ''
        if self.sandwich_mode:
            sandwich = ',sandwich'
        return "Reader(%s,%s,%r%s)" % (
        "[%s]" % (','.join(self.topology)), "[%s]" % (','.join(self.trajectory)), self.window, sandwich)

    def sandwich(self, number=False):
        # generates readers of consecutive sandwich layers
        for nr, trajectory in enumerate(self.trajectory):
            if number:
                yield nr, self.get_single_reader(nr)
            else:
                yield self.get_single_reader(nr)

    def baguette(self, number=False):
        # generates reader of trajectory file(s)
        if number:
            yield 0, self.get_single_reader(0)
        else:
            yield self.get_single_reader(0)

    def strata(self, number=False):
        # generates slices of baquette
        for nr in range(self.number_of_layers()):
            if number:
                yield nr + 1, self.get_single_reader(nr + 1)
            else:
                yield self.get_single_reader(nr + 1)

    def iterate(self, number=False, threads=True):
        # iterates over trajectory readers
        # calls sandwich or baguette
        if self.sandwich_mode:
            return self.sandwich(number=number)
        elif threads and self.threads > 1:
            return self.strata(number=number)
        return self.baguette(number=number)

    def get_single_raw_reader_per_trajectory(self, number):
        # here number is a number of trajectory, starting from 0
        if len(self.topology) == len(self.trajectory):
            return self.engine(self.topology[number],
                               [self.trajectory[number]],
                               number=number,
                               window=Window(0, None, 1))
        else:
            return self.engine(self.topology[0],
                               [self.trajectory[number]],
                               number=number,
                               window=Window(0, None, 1))

    def get_single_reader(self, number):
        # returns single trajectory reader of number
        # is it is already opened it is returned directly
        # if not is opened and then recursive call is executed
        # print "GetSingleReader(%r)" % (number,)
        if self.sandwich_mode:
            if len(self.topology) == 1:
                return self.engine(self.topology[0], self.trajectory[number], number=number, window=self.window)
            else:
                return self.engine(self.topology[number], self.trajectory[number], number=number, window=self.window)
        elif number > 0 and self.threads > 1:
            window = list(self.window.split(self.number_of_layers()))[number - 1]
            return self.engine(self.topology[0], self.trajectory, number=number, window=window)
        else:
            assert number == 0, "Sandwich/baguette mismatch. Try use/not use --sandwich option."
            return self.engine(self.topology[0], self.trajectory, number=0, window=self.window)

    def get_reader_by_id(self, someid):
        # returns trajectory reader of number in someid
        # someid is tuple (number,ix)
        return self.get_single_reader(someid[0])

    def real_number_of_frames(self):
        # returns number of frames in traj files
        # in sandwich it returns number of frames in first layer
        # in baguette it returns total number of frames (so it is also number of frames in layer 0)
        return open_traj_reader(self.get_single_reader(0)).real_number_of_frames()

    def get_edges(self):
        # should return list of edges of trajectory
        # if one trajectory file is used there are two edges begin and end
        # for each additional trajectory file additional edge marking joining point is returned
        # make sense in non sandwich mode only
        # reader has to be reset before and after
        self.reset()
        # to calculate this each individual file should be open as separate trajectory
        # 1. get real number of frames for each part
        rnf = [0]
        for number in range(len(self.trajectory)):
            # if len(rnf)==0:
            #    rnf.append(0)
            # else:
            #    rnf.append(rnf[-1]+1)
            # rnf.append(rnf[-1]-1+open_traj_reader(self.get_single_raw_reader_per_trajectory(number)).real_number_of_frames())
            rnf.append(open_traj_reader(self.get_single_raw_reader_per_trajectory(number)).real_number_of_frames())
        # cumulative sum
        rnf = (np.cumsum(rnf) - 1).tolist()
        # 2. by using window try to calculate where are the edges
        E = []  # list of edges
        for nr, rf in enumerate(self.window.range()):  # iterate over real frames
            if rf >= rnf[0]:
                # this is an edge
                if self.window.step > 1:
                    E.append(nr - 1)
                else:
                    E.append(nr)
                while len(rnf) and rf >= rnf[0]:
                    rnf.pop(0)
        if len(rnf):
            E.append(nr)
        self.reset()
        return E[1:-1]

    def number_of_frames(self, onelayer=False):
        # number of frames in the window times number of layers
        # in baguette number of layers is 1
        if self.sandwich_mode and not onelayer:
            return len(self.trajectory) * self.window.len()
        return self.window.len()

    def number_of_layers(self):
        if self.sandwich_mode:
            return len(self.trajectory)
        nof = self.number_of_frames()
        nol = self.threads * self.threads_multiply
        return nol if nof > nol else max(nof - 1, 1)


# instance of MasterReader
Reader = MasterReader()


def open_traj_reader_engine(ort):
    if ort.engine == 'mda':
        return ReaderTrajViaMDA(ort.topology, ort.trajectory, number=ort.number, window=ort.window)
    raise NotImplementedError


def open_traj_reader(ort):


    if ort.number not in Reader.open_reader_traj:
        Reader.open_reader_traj.update({ort.number: open_traj_reader_engine(ort)})
    return Reader.open_reader_traj[ort.number]


class ReaderAccess(object):
    # ReaderAccess class provides reader property that returns current instance of MasterReader

    def get_reader(self, number):


        if number in Reader.open_reader_traj:
            # print "Getting reader",number,"from opened readers."
            return Reader.open_reader_traj[number]
        Reader.get_single_reader(number).open()
        return self.get_reader(number)

    def get_reader_by_id(self, someid):
        return self.get_reader(someid[0])

    def get_edges(self):


        return Reader.edges


################################################################################
# VdW radii

VdW_radii = {'ac': 2.47,
             'ag': 2.11,
             'al': 1.84,
             'am': 2.44,
             'ar': 1.88,
             'as': 1.85,
             'at': 2.02,
             'au': 2.14,
             'b': 1.92,
             'ba': 2.68,
             'be': 1.53,
             'bi': 2.07,
             'bk': 2.44,
             'br': 1.85,
             'c': 1.7,
             'ca': 2.31,
             'cd': 2.18,
             'ce': 2.42,
             'cf': 2.45,
             'cl': 1.75,
             'cm': 2.45,
             'co': 2.0,
             'cr': 2.06,
             'cs': 3.43,
             'cu': 1.96,
             'dy': 2.31,
             'er': 2.29,
             'es': 2.45,
             'eu': 2.35,
             'f': 1.47,
             'fe': 2.04,
             'fm': 2.45,
             'fr': 3.48,
             'ga': 1.87,
             'gd': 2.34,
             'ge': 2.11,
             'h': 1.1,
             'he': 1.4,
             'hf': 2.23,
             'hg': 2.23,
             'ho': 2.3,
             'i': 1.98,
             'in': 1.93,
             'ir': 2.13,
             'k': 2.75,
             'kr': 2.02,
             'la': 2.43,
             'li': 1.82,
             'lr': 2.46,
             'lu': 2.24,
             'md': 2.46,
             'mg': 1.73,
             'mn': 2.05,
             'mo': 2.17,
             'n': 1.55,
             'na': 2.27,
             'nb': 2.18,
             'nd': 2.39,
             'ne': 1.54,
             'ni': 1.97,
             'no': 2.46,
             'np': 2.39,
             'o': 1.52,
             'os': 2.16,
             'p': 1.8,
             'pa': 2.43,
             'pb': 2.02,
             'pd': 2.1,
             'pm': 2.38,
             'po': 1.97,
             'pr': 2.4,
             'pt': 2.13,
             'pu': 2.43,
             'ra': 2.83,
             'rb': 3.03,
             're': 2.16,
             'rh': 2.1,
             'rn': 2.2,
             'ru': 2.13,
             's': 1.8,
             'sb': 2.06,
             'sc': 2.15,
             'se': 1.9,
             'si': 2.1,
             'sm': 2.36,
             'sn': 2.17,
             'sr': 2.49,
             'ta': 2.22,
             'tb': 2.33,
             'tc': 2.16,
             'te': 2.06,
             'th': 2.45,
             'ti': 2.11,
             'tl': 1.96,
             'tm': 2.27,
             'u': 2.41,
             'v': 2.07,
             'w': 2.18,
             'xe': 2.16,
             'y': 2.32,
             'yb': 2.26,
             'zn': 2.01,
             'zr': 2.23}
'''
Dictionary of VdW radii.

Data taken from L. M. Mentel, mendeleev, 2014. Available at: https://bitbucket.org/lukaszmentel/mendeleev.
Package **mendeleev** is not used because it depends on too many other libraries.
'''


################################################################################
# Reader Traj base class

# class ReaderTraj(ReaderAccess):
class ReaderTraj(object):
    def __init__(self, topology, trajectory,
                 number=None, window=None):
        """
        :param str topology:  Topology file name.
        :param list trajectory: Trajectory file name.
        :param int number: Number of trajectory file.
        :param Window window: Frames window to read.
        :param Reader reader: Parent reader object.

        This is base class for MD data access engines.
        """

        self.topology = topology  # here, topology is only one file
        self.trajectory = trajectory

        if not isinstance(trajectory, list):
            self.trajectory = [t.strip() for t in trajectory.split(pathsep)]

        if not topology or (isinstance(topology,str) and len(topology.strip())==0):
            raise TypeError('Empty topology file')
        # print "ReaderTraj(%r,%r)" % (self.topology,self.trajectory)

        self.number = number
        self.window = window

        self.real_frame_nr = None

        self.trajectory_object = self.open_trajectory()

    def __repr__(self):
        if len(self.trajectory) == 1:
            return "%s(%s,%s,%d,%r)" % (
                self.__class__.__name__, self.topology, self.trajectory[0], self.number, self.window)
        return "%s(%s,%s,%d,%r)" % (
            self.__class__.__name__, self.topology, "[%s]" % (','.join(self.trajectory)), self.number, self.window)

    def iterate_over_frames(self, reverse=False):
        window_len = self.window.len()
        for frame, real_frame in enumerate(self.window.range(reverse=reverse)):
            self.set_real_frame(real_frame)
            if reverse:
                yield window_len - frame - 1
            else:
                yield frame

    def iterate(self, reverse=False):
        return self.iterate_over_frames(reverse=reverse)

    def set_frame(self, frame):
        return self.set_real_frame(self.window.get_real(frame))

    def dump_frames(self, frames, selection=None, filename=None):
        if filename is None:
            filename = create_tmpfile(ext='pdb')
        self.dump_frames_to_file(frames, selection, filename)
        return filename

    def __del__(self):
        self.close_trajectory()

    def number_of_frames(self):
        # should return number of frames in this reader (window)
        return self.window.len()

    def open_trajectory(self):
        # should return any object that can be further used to parse trajectory
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def close_trajectory(self):
        # should close trajectory reader in self.trajectory_object
        # WARNING: This method has to be carefully implemented because it is used by __del__ and
        #          should not emmit any error messages. This is subjet of change.
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def set_real_frame(self, real_frame):
        self.real_frame_nr = real_frame
        # sets real frame
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def real_number_of_frames(self):
        # should return number of frames
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def parse_selection(self, selection):
        # should return selection
        # selection resolution should be atoms
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def select_all(self):
        # should return selection of all atoms
        # selection resolution should be atoms
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def atom_vdw(self, atomid):
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def atom2residue(self, atomid):
        # one atom id to one residue id
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def atoms_positions(self, atomids):
        # atoms ids to coordinates
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def residues_positions(self, resids):
        # residues ids to center of masses coordinates
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def residues_names(self, resids):
        # residues ids to center of masses coordinates
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def atoms_masses(self, atomids):
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def dump_frames_to_file(self, frames, selection, filename):
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")


################################################################################
# ReaderTraj engine MDAnalysis

mda_available_formats = {re.compile('(nc|NC)'): 'nc',
                         re.compile('(prmtop|parmtop|top|PRMTOP|PARMTOP|TOP)'): 'PRMTOP',
                         re.compile('(dcd|DCD)'): 'LAMMPS',
                         re.compile('(psf|PSF)'): 'psf',
                         re.compile('(pdb|PDB)'): 'pdb',
                         re.compile('(crd|CRD)'): 'crd',
                         re.compile('(xtc|XTC)'): 'XTC'}


class ReaderTrajViaMDA(ReaderTraj):

    def open_trajectory(self):
        topology = splitext(self.topology)[1][1:]
        for afk in list(mda_available_formats.keys()):
            if afk.match(topology):
                topology = mda_available_formats[afk]
                break
        trajectory = splitext(self.trajectory[0])[1][1:]
        for afk in list(mda_available_formats.keys()):
            if afk.match(trajectory):
                trajectory = mda_available_formats[afk]
                break
        #print("%"*80)
        #print(topology,trajectory)
        return mda.Universe(self.topology,
                            self.trajectory,
                            topology_format=topology,
                            format=trajectory)

    def close_trajectory(self):
        if hasattr(self, 'trajectory_object'):
            if hasattr(self.trajectory_object, 'trajectory'):
                if hasattr(self.trajectory_object.trajectory, 'close'):
                    self.trajectory_object.trajectory.close()

    def set_real_frame(self, real_frame):
        self.real_frame_nr = real_frame
        self.trajectory_object.trajectory[real_frame]

    def parse_selection(self, selection):
        return AtomSelection({self.number: self.trajectory_object.select_atoms(selection).atoms.ix})

    def select_all(self):
        return self.parse_selection("name *")

    def atom2residue(self, atomid):
        return self.trajectory_object.atoms[atomid].residue.ix

    def atoms_positions(self, atomids):
        # atoms ids to coordinates
        return self.trajectory_object.atoms[atomids].positions

    def residues_positions(self, resids):
        # residues ids to center of masses coordinates
        return (res.atoms.center_of_geometry() for res in
                self.trajectory_object.residues[list(resids)])

    def residues_names(self, resids):
        # residues ids to center of masses coordinates
        for rid in resids:
            yield self.trajectory_object.residues[rid].resname

    def real_number_of_frames(self):
        # should return number of frames
        return len(self.trajectory_object.trajectory)

    def atoms_masses(self, atomids):
        return self.trajectory_object.atoms[atomids].masses

    def atom_vdw(self, atomid):
        element = str(guess_atom_element(self.trajectory_object.atoms[atomid].name)).lower()
        if element in VdW_radii:
            return VdW_radii[element]
        element = str(self.trajectory_object.atoms[atomid].element).lower()
        if element in VdW_radii:
            return VdW_radii[element]
        return 1.4

    def dump_frames_to_file(self, frames, selection, filename):
        mdawriter = mda.Writer(filename, multiframe=True)
        to_dump = self.trajectory_object.select_atoms(selection)
        for frame in frames:
            self.set_frame(frame)
            mdawriter.write(to_dump)
        mdawriter.close()


################################################################################
# Selection objects

class Selection(ReaderAccess):

    def __init__(self, selected):

        self.selected = OrderedDict(selected)
        for number, ids in self.selected.items():
            self.selected[number] = list(map(defaults.int_default, ids))

    def layer(self, number):
        if number in self.selected:
            return self.__class__({number: self.selected[number]})
        return self.__class__({})

    def numbers(self):
        return list(self.selected.keys())

    def ix(self, ix):
        # gets selection of index ix
        ix_current = 0
        for number, ids in self.selected.items():
            if ix_current + len(ids) >= ix + 1:
                # it is here!
                return self.__class__({number: [ids[ix - ix_current]]})  # FIXME: looks like a bug!
            ix_current += len(ids)
        raise IndexError()

    def single(self, selection_id):
        # TODO: suboptimal and probably does not make any sense
        if selection_id[0] in self.selected:
            if selection_id[-1] in self.selected[selection_id[0]]:
                return self.__class__({selection_id[0]: [selection_id[-1]]})
        raise IndexError()

    def len(self):
        _len = 0
        for number, ids in self.selected.items():
            _len += len(ids)
        return _len

    def __len__(self):
        return self.len()

    def add(self, other):

        for number, ids in other.selected.items():
            if number in self.selected:
                self.selected[number] = self.selected[number] + list(ids)
            else:
                self.selected.update({number: ids})

    def remove(self, other):
        """
        Remove all items that exist in other selection.

        :param other: Other selection.
        """
        empty_keys = []
        for number, ids in other.selected.items():
            self.selected[number] = [id_ for id_ in self.selected[number] if id_ not in ids]

            if not self.selected[number]:
                empty_keys.append(number)

        for k in empty_keys:
            del self.selected[k]

    def uniquify(self):

        for number, ids in self.selected.items():
            self.selected[number] = sorted(set(ids))

    def ids(self):
        # report it in this order! always!
        # these are not unique! run run uniqify first!
        for number, ids in self.selected.items():
            for i in ids:
                yield (number, i)

    def coords(self):
        # order of coords should be the same as in ids!
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def center_of_mass(self):
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")


################################################################################

class AtomSelection(Selection):

    def vdw(self):
        for number, ids in self.selected.items():
            for aid in ids:
                yield self.get_reader(number).atom_vdw(aid)

    def residues(self):
        # returns residues selection
        def get_unique_residues():
            for number, ids in self.selected.items():
                number_reader = self.get_reader(number)
                yield number, sorted(set(map(number_reader.atom2residue, ids)))

        return ResidueSelection(get_unique_residues())

    def coords(self):
        for number, ids in self.selected.items():
            number_reader = self.get_reader(number)
            for coord in number_reader.atoms_positions(ids):  # .tolist():
                yield coord

    def center_of_mass(self):
        center = np.zeros(3)
        total_mass = 0
        for number, ids in self.selected.items():
            masses = self.get_reader(number).atoms_masses(ids)
            masses.shape = (len(masses), 1)
            total_mass += sum(masses)
            center += (masses * self.get_reader(number).atoms_positions(ids)).sum(0)
        return center / total_mass

    def contains_residues(self, other_residues, convex_hull=False, convex_hull_inflate=None, map_fun=None,
                          known_true=None):
        # FIXME: known_true slows down!
        # known_true are only ix!
        assert isinstance(other_residues, ResidueSelection)
        if map_fun is None:
            map_fun = map
        # known_true should be list of ids
        if known_true is not None and len(known_true):
            this_ids = list(self.residues().ids())
            kt_list = []
            ch_list = []
            other_coords = []
            for other_id, other_coord in zip(other_residues.ids(), other_residues.coords()):
                if other_id[-1] in known_true:
                    kt_list.append(other_id)
                elif convex_hull:
                    ch_list.append(other_id)
                    other_coords.append(other_coord)
                elif other_id in this_ids:
                    kt_list.append(other_id)
            if convex_hull:
                chull = self.chull(inflate=convex_hull_inflate)
                # ch_results = map_fun(is_point_within_convexhull, izip_longest(other_coords, [], fillvalue=chull))
                ch_results = are_points_within_convexhull(other_coords, chull)
            # final merging loop
            final_results = []
            for other_id in other_residues.ids():
                if other_id in ch_list:
                    final_results.append(ch_results[ch_list.index(other_id)])
                elif other_id in kt_list:
                    final_results.append(True)
                else:
                    final_results.append(False)
            return final_results
        else:
            if convex_hull:
                other_coords = other_residues.coords()
                chull = self.chull(inflate=convex_hull_inflate)
                # return map_fun(is_point_within_convexhull, izip_longest(other_coords, [], fillvalue=chull))
                return are_points_within_convexhull(other_coords, chull)
            else:
                # check if other selection is empty
                this_ids = list(self.residues().ids())
                return [other_id in this_ids for other_id in other_residues.ids()]

    def containing_residues(self, other_residues, *args, **kwargs):
        # Convenience wrapper for contains_residues.
        # def get_res_in_scope(is_res_in_scope, res):
        other_new = OrderedDict()
        for iris, (number, resid) in zip(self.contains_residues(other_residues, *args, **kwargs), other_residues.ids()):
            if iris:
                if number in other_new:
                    other_new[number].append(resid)
                else:
                    other_new.update({number: [resid]})
        return ResidueSelection(other_new)

    def chull(self, inflate=None):
        if self.len() > 3:
            return SciPyConvexHull(list(self.coords()), inflate=inflate)


################################################################################

class ResidueSelection(Selection):

    def coords(self):
        for number, ids in self.selected.items():
            number_reader = self.get_reader(number)
            for coord in number_reader.residues_positions(ids):
                yield coord  # .tolist()

    def names(self):
        for number, ids in self.selected.items():
            for name in self.get_reader(number).residues_names(ids):
                yield name

    def single_residues(self):
        for resid in self.ids():
            yield SingleResidueSelection(resid)


################################################################################

@memory
def coords_range_core(srange, number, rid):
    assert isinstance(srange, SmartRangeIncrement), "Expecting SmartRangeIncrement, got %r" % type(srange)

    @arrayify(shape=(None, 3))
    def coords_range_core_inner(srange, number, rid):


        reader = Reader.get_single_reader(number).open()
        for f in srange.get():
            reader.set_frame(f)
            yield next(reader.residues_positions([rid]))

    return coords_range_core_inner(srange, number, rid)


@arrayify(shape=(None, 3))
def coords_range(srange, number, rid):
    # srange is SmartRangeIncrement, it cannot be anything else
    # wrapper to limit number of calls to coords_range_core
    # CRIC cache is {number:{rid:FramesRangeCollection}}
    logger.debug("CRIC global request %d:%d %s", number, rid, str(srange))
    '''
    if number not in CRIC.cache:
        CRIC.cache.update({number:{}})
        logger.debug("CRIC new number %d",number)
    if rid not in CRIC.cache[number]:
        CRIC.cache[number].update({rid: FramesRangeCollection()})
        logger.debug("CRIC new rid %d",rid)
    '''

    # call get_ranges and stack array, do it in comprehension? nested function?
    def get_coords_from_cache():
        for sr, xr in CRIC.get_frc(number, rid).get_ranges(srange):
            logger.debug("CRIC partial request %d:%d %s", number, rid, str(sr))
            yield coords_range_core(sr, number, rid)[xr, :]

    return np.vstack(tuple(get_coords_from_cache()))


################################################################################

@memory
def smooth_coords_ranges(sranges, number, rid, smooth):
    @tupleify
    def smooth_coords_ranges_inner(sranges, number, rid, smooth):
        # here, sranges are list of SmartRange objects whereas coords_range accepts SmartRangeIncrement
        # first get all coords and make in continous
        def sranges2coords_cont():
            for srange in sranges:
                for srangei in srange.raw:
                    yield coords_range(srangei, number, rid)

        coords_cont = make_default_array(np.vstack([c for c in sranges2coords_cont() if len(c) > 0]))
        # call smooth
        coords_cont = smooth(coords_cont)
        # split coords_cont
        # now lets return tupple of coords
        nr = 0
        for path in sranges:
            if len(path) > 0:
                yield make_default_array(coords_cont[nr:nr + len(path)])
                nr += len(path)
            else:
                # TODO: reuse some kind of empty coords
                yield make_default_array(np.zeros((0, 3)))

    return smooth_coords_ranges_inner(sranges, number, rid, smooth)


################################################################################

class SingleResidueSelection(ReaderAccess):
    def __init__(self, resid):
        super(SingleResidueSelection, self).__init__()
        # where resid is id reported by ResidueSelection and reader is Reader
        # resid is tuple (number,id) number is used to get reader_traj
        self.resid = resid[-1]
        self.number = resid[0]
        # to optionally simplify CRIC store lets read minimal and maximal frame
        # and always ask for full range and then slice it according to frames.
        traj_reader = self.get_reader(self.number)
        # full range always
        self.always_request_frames = SmartRangeIncrement(0,traj_reader.number_of_frames())


    def coords(self, frames):
        if isinstance(frames, SmartRange):
            if len(frames):
                if GCS.cachetype == 'simple':
                    # this is naive version
                    full_coords = coords_range(self.always_request_frames,self.number,self.resid)
                    #return np.vstack([full_coords[list(srange.get()),:] for srange in frames.raw])
                    def get_partial_coords():
                        for sr in frames.raw:
                            lo = sr.first_element() - self.always_request_frames.first_element()
                            hi = lo + len(sr)
                            yield full_coords[lo:hi,:]
                    return np.vstack(list(get_partial_coords()))
                else: # if full
                    # the normal way
                    return np.vstack([coords_range(srange, self.number, self.resid) for srange in frames.raw])
            return self._coords([])
        else:
            return self._coords(frames)

    @arrayify(shape=(None, 3))
    def _coords(self, frames):
        # return coords for frames
        if len(frames):
            # print "Trying to get reader",self.number,"for",len(frames),"frames"
            traj_reader = self.get_reader(self.number)
            for f in frames:
                traj_reader.set_frame(f)
                yield next(traj_reader.residues_positions([self.resid]))

    def coords_smooth(self, sranges, smooth):
        for coord in smooth_coords_ranges(sranges, self.number, self.resid, smooth):
            yield coord


################################################################################

# flag sandwich as imported
GCS.sandwich_import = True
