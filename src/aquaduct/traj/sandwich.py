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


from os.path import splitext
import re
import MDAnalysis as mda
import numpy as np

from collections import OrderedDict
from itertools import izip_longest


from aquaduct.utils.helpers import is_iterable
from aquaduct.geom.convexhull import SciPyConvexHull,is_point_within_convexhull


class Window(object):
    def __init__(self,start,stop,step):
        self.start = start
        self.stop = stop
        self.step = step

    def __repr__(self):
        return "Window(%r:%r:%r)" % (self.start,self.stop,self.step)

    def range(self):
        return xrange(self.start,self.stop,self.step)

    def get_real(self,frame):
        # recalculates real frame
        return self.start + frame*self.step

    def len(self):
        return (self.start - self.stop) / self.step

# engine problem
# Two options are currently available:
# MDAnalysis
# MDTraj


class Reader(object):

    def __init__(self,topology,trajectory,
                 window=None,
                 sandwich=False):
        '''
        :param str topology:  Topology file name.
        :param list trajectory: List of trajectories. Each element is a fine name.
        :param Window window: Frames window to read.
        :param bool sandwich: Flag for setting sandwitch mode.

        Window is propagated to ReaderTraj.
        '''

        self.topology = topology
        if not isinstance(trajectory,list):
            trajectory = [trajectory]
        self.trajectory = trajectory

        self.window = window
        self.sandwich_mode = sandwich

        # mda
        self.engine = ReaderTrajViaMDA

        self.open_reader_traj = {}
        self.correct_window()
        self.open_reader_traj = {} # clear that

    def correct_window(self):
        # correct window!
        N = self.real_number_of_frames()
        start,stop,step = self.window.start,self.window.stop,self.window.step
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
        self.window = Window(start,stop,step)

    def __repr__(self):
        sandwich = ''
        if self.sandwich_mode:
            sandwich = ',sandwich'
        return "Reader(%s,%s,%r%s)" % (self.topology,"[%s]" % (','.join(self.trajectory)),self.window,sandwich)

    def sandwich(self):
        for nr,trajectory in enumerate(self.trajectory):
            yield self.get_single_reader(nr)

    def baguette(self):
        yield self.get_single_reader(0)

    def iterate(self):
        if self.sandwich_mode:
            return self.sandwich()
        return self.baguette()

    def get_single_reader(self,number):
        if self.open_reader_traj.has_key(number):
            return self.open_reader_traj[number]
        else:
            if self.sandwich_mode:
                self.open_reader_traj.update({number:self.engine(self.topology,self.trajectory[number],number=number,window=self.window,reader=self)})
            else:
                assert number == 0
                self.open_reader_traj.update({0: self.engine(self.topology,self.trajectory[number],number=0,window=self.window,reader=self)})
        return self.get_single_reader(number)

    def real_number_of_frames(self):
        return self.get_single_reader(0).real_number_of_frames()

    def number_of_frames(self):
        if self.sandwich():
            return len(self.trajectory)*self.window.len()
        return self.window.len()

class ReaderTraj(object):
    def __init__(self,topology,trajectory,
                 number=None,window=None,reader=None):
        '''
        :param str topology:  Topology file name.
        :param list trajectory: Trajectory file name.
        :param int number: Number of trajectory file.
        :param Window window: Frames window to read.
        :param Reader reader: Parent reader object.
        '''

        self.topology = topology
        self.trajectory = trajectory
        if not isinstance(trajectory,list):
            self.trajectory = [trajectory]

        self.number = number
        self.window = window
        assert isinstance(reader,Reader)
        self.reader = reader

        self.real_frame_nr = None

        self.trajectory_object = self.open_trajectory()

    def __repr__(self):
        if len(self.trajectory) == 1:
            return "%s(%s,%s,%d,%r)" % (self.__class__.__name__,self.topology,self.trajectory[0],self.number,self.window)
        return "%s(%s,%s,%d,%r)" % (self.__class__.__name__,self.topology,"[%s]" % (','.join(self.trajectory)),self.number,self.window)

    def iterate_over_frames(self):
        for frame,real_frame in enumerate(self.window.range()):
            self.set_real_frame(real_frame)
            yield frame

    def set_frame(self,frame):
        return self.set_real_frame(self.window.get_real(frame))

    def open_trajectory(self):
        # should return any object that can be further used to parse trajectory
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def set_real_frame(self,real_frame):
        self.real_frame_nr = real_frame
        # sets real frame
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def parse_selection(self, selection):
        # should return selection
        # selection resolution should be atoms
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")


    def atom2residue(self,atomid):
        # one atom id to one residue id
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def atoms_positions(self,atomids):
        # atoms ids to coordinates
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def residues_positions(self,resids):
        # residues ids to center of masses coordinates
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def atoms_masses(self,atomids):
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")


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

    def set_real_frame(self,real_frame):
        self.real_frame_nr = real_frame
        self.trajectory_object.trajectory[real_frame]


    def parse_selection(self, selection):
        return AtomSelection({self.number:self.trajectory_object.select_atoms(selection).atoms.ix},
                              reader=self.reader)

    def atom2residue(self,atomid):
        return self.trajectory_object.atoms[atomid].residue.ix

    def atoms_positions(self,atomids):
        # atoms ids to coordinates
        return self.trajectory_object.atoms[atomids].positions

    def residues_positions(self,resids):
        # residues ids to center of masses coordinates
        for rid in resids:
            yield self.trajectory_object.residues[[rid]].center_of_mass()

    def real_number_of_frames(self):
        # should return number of frames
        return len(self.trajectory_object.trajectory)

    def atoms_masses(self,atomids):
        return self.trajectory_object.atoms[atomids].masses


class Selection(object):

    def __init__(self,selected,reader=None):

        assert isinstance(reader,Reader)

        self.selected = OrderedDict(selected)
        for number, ids in self.selected.iteritems():
            self.selected[number] = list(ids)

        self.reader = reader


    def get_reader(self,number):
        return self.reader.get_single_reader(number)

    def add(self,other):

        for number,ids in other.selected.iteritems():
            if self.selected.has_key(number):
                self.selected[number] = self.selected[number]+list(ids)
            else:
                self.selected.update({number:ids})

    def uniqify(self):

        for number, ids in self.selected.iteritems():
            self.selected[number] = sorted(set(ids))

    def ids(self):
        # reportet in this order! always!
        # these are not unique! run run uniqify first!
        for number, ids in self.selected.iteritems():
            for i in ids:
                yield (number,i)

    def coords(self):
        # order of coords should be the same as in ids!
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def center_of_mass(self):
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

class AtomSelection(Selection):

    def residues(self):
        # returns residues selection
        def get_unique_residues():
            for number, ids in self.selected.iteritems():
                yield number,sorted(set(map(self.get_reader(number).atom2residue,ids)))
        return ResidueSelection(get_unique_residues(),reader=self.reader)

    def coords(self):
        for number, ids in self.selected.iteritems():
            for coord in self.get_reader(number).atoms_positions(ids):
                yield coord.tolist()

    def center_of_mass(self):
        center = np.zeros(3)
        total_mass = 0
        for number, ids in self.selected.iteritems():
            masses = self.get_reader(number).atoms_masses(ids)
            masses.shape = (len(masses),1)
            total_mass += sum(masses)
            center += (masses*self.get_reader(number).atoms_positions(ids)).sum(0)
        return center/total_mass

    def contains_residues(self,other_residues,convex_hull=False,map_fun=None,known_true=None):
        assert isinstance(other_residues,ResidueSelection)
        if map_fun is None:
            map_fun = map
        # known_true should be list of ids
        if known_true is not None and len(known_true):
            this_ids = self.residues().ids
            kt_list = []
            ch_list = []
            for other_id in other_residues.ids():
                if other_id in known_true:
                    kt_list.append(other_id)
                elif convex_hull:
                    ch_list.append(other_id)
                elif other_id in this_ids:
                    kt_list.append(other_id)
            if convex_hull:
                other_coords = other_residues.coords()
                chull =  SciPyConvexHull(self.coords())
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
                other_coords = list(other.center_of_mass_of_residues())
                chull = self.get_convexhull_of_atom_positions()
                return map_fun(is_point_within_convexhull, izip_longest(other_coords, [], fillvalue=chull))
            else:
                # check if other selection is empty
                if other.unique_resids_number() == 0:
                    return []
                this_uids = self.unique_resids(ikwid=True)
                return [res_other.unique_resids(ikwid=True) in this_uids for res_other in other.iterate_over_residues()]



class ResidueSelection(Selection):

    def coords(self):
        for number, ids in self.selected.iteritems():
            for coord in self.get_reader(number).residues_positions(ids):
                yield coord.tolist()


class SingleResidueSelection(object):
    def __init__(self,resid,reader):
        # where resid is id reported by ResidueSelection and reader is Reader
        # resid is tuple (number,id) number is used to get reader_traj
        self.resid = resid[-1]
        self.number = resid[0]
        self.reader_traj = reader.get_single_reader(self.number)

    def coords(self,frames):
        # return coords for frames
        for f in frames:
            self.reader_traj.set_frame(f)
            yield self.reader_traj.residues_positions([self.resid]).next()
