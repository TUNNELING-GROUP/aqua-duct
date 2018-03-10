# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2018  Tomasz Magdziarz <info@aquaduct.pl>
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


import re

from os.path import splitext
from collections import OrderedDict
from itertools import izip_longest, imap

import numpy as np

import MDAnalysis as mda
from MDAnalysis.topology.core import guess_atom_element

from aquaduct.utils.helpers import is_iterable
from aquaduct.geom.convexhull import SciPyConvexHull, is_point_within_convexhull
from aquaduct.utils.helpers import arrayify, SmartRange, create_tmpfile, \
                                   tupleify, SmartRangeIncrement
from aquaduct.utils.maths import make_default_array
from aquaduct.apps.data import GCS,CRIC
from aquaduct.utils.maths import defaults
from aquaduct import logger

################################################################################
# memory decorator

if GCS.cachedir:
    from joblib import Memory

    memory_cache = Memory(cachedir=GCS.cachedir, verbose=0) # mmap have to be switched off, otherwise smoothing does not work properly
    #memory_cache = Memory(cachedir=GCS.cachedir, mmap_mode='r', verbose=0)
    memory = memory_cache.cache
elif GCS.cachemem:
    from functools import wraps


    # TODO: rework this to adapt it accordingly, moce it to utils?
    # https://medium.com/@nkhaja/memoization-and-decorators-with-python-32f607439f84
    def memory(func):
        cache = func.cache = {}

        @wraps(func)
        def memoized_func(*args, **kwargs):
            key = ','.join(map(str,args)) + '&' + ','.join(map(lambda kv: ':'.join(map(str,kv)),kwargs.iteritems()))
            logger.debug('Looking for cache key %s' % key)
            if key not in cache:
                cache[key] = func(*args, **kwargs)
                logger.debug("New key added to cache.")
            return cache[key]

        return memoized_func
else:
    from aquaduct.utils.helpers import noaction as memory


################################################################################
# trajectory window object

class Window(object):
    def __init__(self, start, stop, step):
        self.start = start
        self.stop = stop
        self.step = step

    def __repr__(self):
        return "Window(%r:%r:%r)" % (self.start, self.stop, self.step)

    def range(self):
        # returns range object
        return xrange(self.start, self.stop + 1, self.step)

    def get_real(self, frame):
        # recalculates real frame
        return self.start + frame * self.step

    def len(self):
        # lenght of window
        return len(self.range())


################################################################################
# MasterReader

# engine problem
# Two options are currently available:
# MDAnalysis
# MDTraj

class MasterReader(object):
    # only one MasterReader object can be used
    # it does not use ANY directo call to ANY MD access software

    open_reader_traj = {}

    topology = ''
    trajectory = ['']
    window = None

    sandwich_mode = None
    engine_name = 'mda'

    def __call__(self, topology, trajectory, window=None, sandwich=False):
        '''
        :param str topology:  Topology file name.
        :param list trajectory: List of trajectories. Each element is a fine name.
        :param Window window: Frames window to read.
        :param bool sandwich: Flag for setting sandwitch mode.
        '''

        self.topology = topology
        if not isinstance(trajectory, list):
            trajectory = [trajectory]
        self.trajectory = trajectory

        self.window = window
        self.sandwich_mode = sandwich

        self.correct_window()  # this corrects window
        self.open_reader_traj = {}  # aster window correction clear all opened trajs

    def __getstate__(self):
        # if pickle dump is required, this will not be used in the future
        # do not pass open_reader_traj
        return dict(((k, v) for k, v in self.__dict__.iteritems() if k not in ['open_reader_traj']))

    def __setstate__(self, state):
        # if pickle dump is required, this will not be used in the future
        self.__dict__ = state
        self.open_reader_traj = {}

    @property
    def engine(self):
        # returns engine used for accessing trajectory
        # this returns object that inherits from ReaderTraj
        if self.engine_name == 'mda':
            return ReaderTrajViaMDA
        raise NotImplementedError

    def correct_window(self):
        # corrects window object
        # should be called at the initialization stage which is executed by calling the object
        # side effect is that it opens at least one trajectory file
        N = self.real_number_of_frames()
        start, stop, step = self.window.start, self.window.stop, self.window.step
        if start is None:
            start = 0
        if stop is None:
            stop = N - 1
        if step is None:
            step = 1
        if start < 0:
            start = 0
        if stop < 0:
            stop = 0
        if start > N:
            start = N - 1
        if stop > N:
            stop = N - 1
        self.window = Window(start, stop, step)

    def __repr__(self):
        sandwich = ''
        if self.sandwich_mode:
            sandwich = ',sandwich'
        return "Reader(%s,%s,%r%s)" % (self.topology, "[%s]" % (','.join(self.trajectory)), self.window, sandwich)

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

    def iterate(self, number=False):
        # iterates over trajectory readers
        # calls sandwich or baguette
        if self.sandwich_mode:
            return self.sandwich(number=number)
        return self.baguette(number=number)

    def get_single_reader(self, number):
        # returns single trajectory reader of number
        # is it is already opened it is returned directly
        # if not is opened and then recursive call is executed
        if self.open_reader_traj.has_key(number):
            return self.open_reader_traj[number]
        else:
            if self.sandwich_mode:
                self.open_reader_traj.update(
                    {number: self.engine(self.topology, self.trajectory[number], number=number, window=self.window)})
            else:
                assert number == 0
                self.open_reader_traj.update(
                    {0: self.engine(self.topology, self.trajectory, number=0, window=self.window)})
        return self.get_single_reader(number)

    def get_reader_by_id(self, someid):
        # returns trajectory reader of number in someid
        # someid is tuple (number,ix)
        return self.get_single_reader(someid[0])

    def real_number_of_frames(self):
        # returns number of frames in traj files
        # in sandwich it returns number of frames in first layer
        # in baguette it returns total number of frames (so it is also number of frames in layer 0)
        return self.get_single_reader(0).real_number_of_frames()

    def number_of_frames(self, onelayer=False):
        # number of frames in the window times number of layers
        # in baguette number of layers is 1
        if self.sandwich_mode and not onelayer:
            return len(self.trajectory) * self.window.len()
        return self.window.len()

    def number_of_layers(self):
        if self.sandwich_mode:
            return len(self.trajectory)
        return 1

# instance of MasterReader
Reader = MasterReader()


class ReaderAccess(object):
    # ReaderAccess class provides reader property that returns current instance of MasterReader

    @property
    def reader(self):
        return Reader


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

class ReaderTraj(ReaderAccess):
    def __init__(self, topology, trajectory,
                 number=None, window=None):
        '''
        :param str topology:  Topology file name.
        :param list trajectory: Trajectory file name.
        :param int number: Number of trajectory file.
        :param Window window: Frames window to read.
        :param Reader reader: Parent reader object.

        This is base class for MD data access engines.
        '''

        self.topology = topology
        self.trajectory = trajectory
        if not isinstance(trajectory, list):
            self.trajectory = [trajectory]

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

    def iterate_over_frames(self):
        for frame, real_frame in enumerate(self.window.range()):
            self.set_real_frame(real_frame)
            yield frame

    def set_frame(self, frame):
        return self.set_real_frame(self.window.get_real(frame))

    def dump_frames(self, frames, selection=None, filename=None):
        if filename is None:
            filename = create_tmpfile(ext='pdb')
        self.dump_frames_to_file(frames, selection, filename)
        return filename

    def __del__(self):
        self.close_trajectory()

    def open_trajectory(self):
        # should return any object that can be further used to parse trajectory
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def close_trajectory(self):
        # should close trajectory reader in self.trajectory_object
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
                         re.compile('(parmtop|top|PARMTOP|TOP)'): 'parmtop',
                         re.compile('(dcd|DCD)'): 'LAMMPS',
                         re.compile('(psf|PSF)'): 'psf',
                         re.compile('(pdb|PDB)'): 'pdb',
                         re.compile('(crd|CRD)'): 'crd',
                         re.compile('(xtc|XTC)'): 'XTC'}


class ReaderTrajViaMDA(ReaderTraj):

    def open_trajectory(self):
        topology = splitext(self.topology)[1][1:]
        for afk in mda_available_formats.keys():
            if afk.match(topology):
                topology = mda_available_formats[afk]
                break
        trajectory = splitext(self.trajectory[0])[1][1:]
        for afk in mda_available_formats.keys():
            if afk.match(trajectory):
                trajectory = mda_available_formats[afk]
                break
        return mda.Universe(self.topology,
                            self.trajectory,
                            topology_format=topology,
                            format=trajectory)

    def close_trajectory(self):
        self.trajectory_object.trajectory.close()

    def set_real_frame(self, real_frame):
        self.real_frame_nr = real_frame
        self.trajectory_object.trajectory[real_frame]

    def parse_selection(self, selection):
        return AtomSelection({self.number: self.trajectory_object.select_atoms(selection).atoms.ix})

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
        for number, ids in self.selected.iteritems():
            self.selected[number] = list(imap(defaults.int_default,ids))

    def layer(self, number):
        if self.selected.has_key(number):
            return self.__class__({number: self.selected[number]})
        return self.__class__({})

    def numbers(self):
        return self.selected.keys()

    def ix(self, ix):
        # gets selection of index ix
        ix_current = 0
        for number, ids in self.selected.iteritems():
            if ix_current + len(ids) >= ix + 1:
                # it is here!
                return self.__class__({number: [ids[ix - ix_current]]})
            ix_current += len(ids)
        raise IndexError()

    def len(self):
        _len = 0
        for number, ids in self.selected.iteritems():
            _len += len(ids)
        return _len

    def get_reader(self, number):
        return self.reader.get_single_reader(number)

    def add(self, other):

        for number, ids in other.selected.iteritems():
            if self.selected.has_key(number):
                self.selected[number] = self.selected[number] + list(ids)
            else:
                self.selected.update({number: ids})

    def uniquify(self):

        for number, ids in self.selected.iteritems():
            self.selected[number] = sorted(set(ids))

    def ids(self):
        # report it in this order! always!
        # these are not unique! run run uniqify first!
        for number, ids in self.selected.iteritems():
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
        for number, ids in self.selected.iteritems():
            for aid in ids:
                yield self.get_reader(number).atom_vdw(aid)

    def residues(self):
        # returns residues selection
        def get_unique_residues():
            for number, ids in self.selected.iteritems():
                number_reader = self.get_reader(number)
                yield number, sorted(set(map(number_reader.atom2residue, ids)))

        return ResidueSelection(get_unique_residues())

    def coords(self):
        for number, ids in self.selected.iteritems():
            number_reader = self.get_reader(number)
            for coord in number_reader.atoms_positions(ids):#.tolist():
                yield coord

    def center_of_mass(self):
        center = np.zeros(3)
        total_mass = 0
        for number, ids in self.selected.iteritems():
            masses = self.get_reader(number).atoms_masses(ids)
            masses.shape = (len(masses), 1)
            total_mass += sum(masses)
            center += (masses * self.get_reader(number).atoms_positions(ids)).sum(0)
        return center / total_mass

    def contains_residues(self, other_residues, convex_hull=False, map_fun=None, known_true=None):
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
                chull = self.chull()
                ch_results = map_fun(is_point_within_convexhull, izip_longest(other_coords, [], fillvalue=chull))
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
                chull = self.chull()
                return map_fun(is_point_within_convexhull, izip_longest(other_coords, [], fillvalue=chull))
            else:
                # check if other selection is empty
                this_ids = list(self.residues().ids())
                return [other_id in this_ids for other_id in other_residues.ids()]

    def containing_residues(self, other_residues, *args, **kwargs):
        # Convienience wrapper for contains_residues.
        # def get_res_in_scope(is_res_in_scope, res):
        other_new = OrderedDict()
        for iris, (number, resid) in zip(self.contains_residues(other_residues, *args, **kwargs), other_residues.ids()):
            if iris:
                if other_new.has_key(number):
                    other_new[number].append(resid)
                else:
                    other_new.update({number: [resid]})
        return ResidueSelection(other_new)

    def chull(self):
        return SciPyConvexHull(list(self.coords()))


################################################################################

class ResidueSelection(Selection):

    def coords(self):
        for number, ids in self.selected.iteritems():
            number_reader = self.get_reader(number)
            for coord in number_reader.residues_positions(ids):
                yield coord #.tolist()

    def names(self):
        for number, ids in self.selected.iteritems():
            for name in self.get_reader(number).residues_names(ids):
                yield name

    def single_residues(self):
        for resid in self.ids():
            yield SingleResidueSelection(resid)


################################################################################

@memory
@arrayify(shape=(None, 3))
def coords_range_core(srange, number, rid):
    reader = Reader.get_single_reader(number)
    for f in srange.get():
        reader.set_frame(f)
        yield reader.residues_positions([rid]).next()



def coords_range(srange, number, rid):
    # srange is SmartRangeIncrement, it cannot be anything else
    # wrapper to limit number of calls to coords_range_core
    # CRIC cache is {number:{rid:FramesRangeCollection}}
    logger.debug("CRIC request %d:%d %s",number,rid,str(srange))
    if number not in CRIC.cache:
        CRIC.cache.update({number:{}})
        logger.debug("CRIC new number %d",number)
    if rid not in CRIC.cache[number]:
        CRIC.cache[number].update({rid:FramesRangeCollection()})
        logger.debug("CRIC new rid %d",rid)
    # call get_ranges and stack array, do it in comprehension? nested function?
    def get_coords_from_cache():
        for sr,xr in CRIC.cache[number][rid].get_ranges(srange):
            yield coords_range_core(sr,number,rid)[xr,:]
    return np.vstack(get_coords_from_cache())


class FramesRangeCollection(object):
    # currently it is assumed that samrt ranges increments only are possible
    def __init__(self):
        self.collection = [] # order on this list does matter!

    def append(self,srange):
        if not len(self.collection):
            self.collection.append(srange)
            logger.debug("FRC append first srange %s",str(srange))
            return
        # there are V cases:
        #           |------|            sr
        # 1 |---|                       it is before sr
        # 2      |-----|                it overlaps with sr but begins before
        # 3          |----|             it is contained in sr
        # 4              |----|         it overlaps with sr but ends after
        # 24     |------------|         it overlabs with sr in 2 and 4 way
        # 5                  |----|     it is after sr
        while (srange is not None):
            for nr,sr in enumerate(self.collection):
                # sr
                if sr.overlaps_mutual(srange):# or srange.overlaps(sr):
                    if sr.contains(srange):
                        # case 3
                        srange = None
                        break
                    if srange.first_element() < sr.first_element():
                        # case 2
                        self.collection.insert(nr,SmartRangeIncrement(srange.first_element(),sr.first_element()-srange.first_element()))
                        logger.debug("FRC case 2 insert srange %s at %d",str(SmartRangeIncrement(srange.first_element(),sr.first_element()-srange.first_element())),nr)
                        if srange.last_element() > sr.last_element():
                            # case 24
                            srange = SmartRangeIncrement(sr.last_element()+1,srange.last_element()-sr.last_element())
                            break
                        else:
                            srange = None
                        break
                    if srange.last_element() > sr.last_element():
                        # case 4
                        srange = SmartRangeIncrement(sr.last_element()+1,srange.last_element()-sr.last_element())
                        continue
                else:
                    if srange.last_element() < sr.first_element():
                        # case 1: insert it before sr
                        self.collection.insert(nr,srange)
                        logger.debug("FRC case 1 insert srange %s at %d",str(srange),nr)
                        srange = None
                        break
                    if srange.first_element() > sr.last_element():
                        # case 5: do nothing
                        continue
            # if something is left append it to the end
            if srange is not None and nr == len(self.collection) - 1:
                self.collection.append(srange)
                logger.debug("FRC append remaining srage %s",str(srange))
                srange = None

    def get_ranges(self,srange):
        # yield sranges from collection and appropriate ranges for these sranges
        # assumes append was already called? call it!
        # after it is called only case 3 or 4 is possible or no overlap at all
        self.append(srange)
        for sr in self.collection:
            if sr.overlaps(srange):
                if sr.contains(srange):
                    # case 3
                    yield sr,xrange(srange.first_element()-sr.first_element(),srange.first_element()-sr.first_element()+len(srange))
                    srange = None
                    break
                # case 4
                yield sr,xrange(srange.first_element()-sr.first_element(),srange.first_element()-sr.first_element()+sr.last_element()-srange.first_element()+1)
                srange = SmartRangeIncrement(sr.last_element()+1,srange.last_element()-sr.last_element())






################################################################################

@memory
@tupleify
def smooth_coords_ranges(sranges,number,rid,smooth):
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


################################################################################

class SingleResidueSelection(ReaderAccess):
    def __init__(self, resid):
        super(SingleResidueSelection, self).__init__()
        # where resid is id reported by ResidueSelection and reader is Reader
        # resid is tuple (number,id) number is used to get reader_traj
        self.resid = resid[-1]
        self.number = resid[0]

    def get_reader(self):
        return self.reader.get_single_reader(self.number)

    def coords(self, frames):
        if isinstance(frames, SmartRange):
            if len(frames):
                return np.vstack([coords_range(srange, self.number, self.resid) for srange in frames.raw])
            return self._coords([])
        else:
            return self._coords(frames)

    @arrayify(shape=(None, 3))
    def _coords(self, frames):
        # return coords for frames
        for f in frames:
            self.get_reader().set_frame(f)
            yield self.get_reader().residues_positions([self.resid]).next()

    def coords_smooth(self,sranges,smooth):
        for coord in smooth_coords_ranges(sranges, self.number, self.resid, smooth):
            yield coord

################################################################################
