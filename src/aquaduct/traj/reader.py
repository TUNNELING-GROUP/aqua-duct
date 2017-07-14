# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2016-2017  Tomasz Magdziarz, Alicja Płuciennik, Michał Stolarczyk <info@aquaduct.pl>
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
from MDAnalysis.topology.core import guess_atom_element

from aquaduct.traj.selections import SelectionMDA

mda_available_formats = {re.compile('(nc|NC)'): 'nc',
                         re.compile('(parmtop|top|PARMTOP|TOP)'): 'parmtop',
                         re.compile('(dcd|DCD)'): 'LAMMPS',
                         re.compile('(psf|PSF)'): 'psf',
                         re.compile('(pdb|PDB)'): 'pdb',
                         re.compile('(xtc|XTC)'): 'XTC'}


class Reader(object):
    def __init__(self, topology, trajectory, window=None):
        assert isinstance(topology, str)
        if not isinstance(trajectory, str):
            for trj in trajectory:
                assert isinstance(trj, str)

        self.topology_file_name = topology
        self.trajectory_file_name = trajectory

        self.trajectory_object = self.open_trajectory()

        # if window is set then behaviour of iteration over frames related methods is different
        self.frames_window = window
        self.frames_window_list = None

        # check window
        lo = min(self.get_start_frame(),self.get_stop_frame())
        hi = max(self.get_start_frame(), self.get_stop_frame())
        assert lo >=0 , 'Wrong frames window definition, wrong low bound %d; should be >= 0' % lo
        assert hi <= self.real_number_of_frames, 'Wrong frames window definition, wrong high bound %d; should be <= %d' % (hi,self.real_number_of_frames-1)

        self.set_real_frame(0)  # ensure that by default it starts from frame 0


    def open_trajectory(self):
        # should return any object that can be further used to parse trajectory
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    @property
    def real_number_of_frames(self):
        # should return real number of frames
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def next_frame(self):
        # should set next frame or raise StopIteration
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def set_real_frame(self, frame):
        # should set current frame according to window
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    ###########################################################################
    # window params

    def get_start_frame(self):
        if self.frames_window is None:
            return 0
        if self.frames_window.start is None:
            return 0
        return self.frames_window.start

    def get_step_frame(self):
        if self.frames_window is None:
            return 1
        if self.frames_window.step is None:
            return 1
        return self.frames_window.step

    def get_stop_frame(self):
        if self.frames_window is None:
            return self.real_number_of_frames - 1
        if self.frames_window.stop is None:
            return self.real_number_of_frames - 1
        return self.frames_window.stop

    def get_window_frame_range(self):
        return xrange(self.get_start_frame(),self.get_stop_frame()+1,self.get_step_frame())
        #return xrange(self.get_start_frame(),self.get_stop_frame(),self.get_step_frame())

    ###########################################################################
    # window dependent

    @property
    def number_of_frames(self):
        if self.frames_window is None:
            return self.real_number_of_frames
        start = self.get_start_frame()
        stop = self.get_stop_frame()
        step = self.get_step_frame()
        return (abs(stop-start)+1)/step

    def iterate_over_frames(self):
        # should return list of frames ids or generator returning such a list, and should set appropriate frame
        # appropriate frame means current frame that can be recalculated to real frame (window)
        for current_frame,real_frame in enumerate(self.get_window_frame_range()):
            self.set_real_frame(real_frame)
            yield current_frame

    def cf2rf(self,cf):
        # current frame to real frame
        if self.frames_window_list is None:
            self.frames_window_list = list(self.get_window_frame_range())
        return self.frames_window_list[cf]

    def set_current_frame(self, frame):
        # frame is current frame, derrive real frame
        return self.set_real_frame(self.cf2rf(frame))

    ###########################################################################

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

    ############################################################################
    # window dependent


    ############################################################################

    def set_real_frame(self, frame):
        self.trajectory_object.trajectory[frame]

    @property
    def real_number_of_frames(self):
        # should return number of frames
        return len(self.trajectory_object.trajectory)

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

def atom2vdw_radius(atom):
    '''
    Function tries to guess atom element and checks if it is in :attr:`VdW_radii` dictionary. If it fails 1.4 is returned.
    Guessing is done twice:

    #. Function :func:`MDAnalysis.topology.core.guess_atom_element` is used.
    #. :attr:`MDAnalysis.core.AtomGroup.Atom.element` is used.

    :param MDAnalysis.core.AtomGroup.Atom atom: Atom of interest.
    :rtype: float
    :return: VdW radius.
    '''
    element = str(guess_atom_element(atom.name)).lower()
    if element in VdW_radii:
        return VdW_radii[element]
    element = str(atom.element).lower()
    if element in VdW_radii:
        return VdW_radii[element]
    return 1.4
