# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2016  Tomasz Magdziarz, Alicja Płuciennik, Michał Stolarczyk <info@aquaduct.pl>
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

import MDAnalysis as mda

from aquaduct.traj.selections import SelectionMDA

mda_available_formats = {re.compile('(nc|NC)'): 'nc',
                         re.compile('(parmtop|top|PARMTOP|TOP)'): 'parmtop',
                         re.compile('(dcd|DCD)'): 'LAMMPS',
                         re.compile('(psf|PSF)'): 'psf',
                         re.compile('(pdb|PDB)'): 'pdb'}


class Reader(object):
    def __init__(self, topology, trajectory):
        assert isinstance(topology, str)
        if not isinstance(trajectory, str):
            for trj in trajectory:
                assert isinstance(trj, str)

        self.topology_file_name = topology
        self.trajectory_file_name = trajectory

        self.trajectory_object = self.open_trajectory()
        self.set_current_frame(0)  # ensure that by default it starts from frame 0

    def open_trajectory(self):
        # should return any object that can be further used to parse trajectory
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    @property
    def number_of_frames(self):
        # should return number of frames
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def set_current_frame(self, frame):
        # should set current frame
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def next_frame(self):
        # should set next frame or raise StopIteration
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def iterate_over_frames(self):
        # should return list of frames ids or generator returning such a list, and should set appropriate frame
        current_frame = 0
        try:
            self.set_current_frame(current_frame)
            while True:
                yield current_frame
                current_frame += 1
                self.next_frame()
        except StopIteration:
            pass

    def parse_selection(self, selection):
        # should parse and return selection object
        # is this object updated automatically depends on the particular implementation of Reader class
        # in particular MDA updates postions of atoms accoriding to current frame, to renew selection one has to parse it again
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def select_resnum(self, resnum):
        # should return selection object
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def select_multiple_resnum(self, resnum_list):
        # should return selection object
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")


class ReadViaMDA(Reader):
    @property
    def number_of_frames(self):
        # should return number of frames
        return len(self.trajectory_object.trajectory)

    def set_current_frame(self, frame):
        self.trajectory_object.trajectory[frame]

    def next_frame(self):
        try:
            self.trajectory_object.trajectory.next()
        except:
            raise StopIteration

    def parse_selection(self, selection):
        return SelectionMDA(self.trajectory_object.select_atoms(selection))

    def select_resnum(self, resnum):
        assert isinstance(resnum, (int, long))
        return self.parse_selection("resnum %d" % resnum)

    def select_multiple_resnum(self, resnum_list):
        assert isinstance(resnum_list, list)
        selection = None
        for resnum in resnum_list:
            if selection is None:
                selection = self.select_resnum(resnum)
            else:
                selection += self.select_resnum(resnum)
        return SelectionMDA(selection)

    def __enter__(self):
        return self

    def __exit__(self, typ, value, traceback):
        self.trajectory_object.trajectory.close()

    def open_trajectory(self):
        topology = splitext(self.topology_file_name)[1][1:]
        for afk in mda_available_formats.keys():
            if afk.match(topology):
                topology = mda_available_formats[afk]
                break
        trajectory = splitext(self.trajectory_file_name[0])[1][1:]
        for afk in mda_available_formats.keys():
            if afk.match(trajectory):
                trajectory = mda_available_formats[afk]
                break
        return mda.Universe(self.topology_file_name,
                            self.trajectory_file_name,
                            topology_format=topology,
                            format=trajectory)


class ReadAmberNetCDFviaMDA(ReadViaMDA):
    def open_trajectory(self):
        return mda.Universe(self.topology_file_name,
                            self.trajectory_file_name,
                            topology_format="prmtop", format="nc")


class ReadDCDviaMDA(ReadViaMDA):
    def open_trajectory(self):
        return mda.Universe(self.topology_file_name,
                            self.trajectory_file_name,
                            topology_format="psf", format="dcd")
