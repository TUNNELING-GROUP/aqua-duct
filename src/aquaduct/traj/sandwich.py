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


class Window(object):
    def __init__(self,start,stop,step):
        self.start = start
        self.stop = stop
        self.step = step

    def __repr__(self):
        return "Window(%r:%r:%r)" % (self.start,self.stop,self.step)

    def range(self):
        return xrange(self.start,self.stop,self.step)

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
        self.trajectory = trajectory

        self.window = window
        self.sandwich_mode = sandwich

        # mda
        self.engine = ReaderTrajViaMDA

    def __repr__(self):
        sandwich = ''
        if self.sandwich_mode:
            sandwich = ',sandwich'
        return "Reader(%s,%s,%r%s)" % (self.topology,"[%s]" % (','.join(self.trajectory)),self.window,sandwich)

    def sandwich(self):
        for nr,trajectory in enumerate(self.trajectory):
            yield self.engine(self.topology,trajectory,
                              number=nr,
                              window=self.window)

    def baguette(self):
        yield self.engine(self.topology,self.trajectory,
                          number=0,
                          window=self.window)

    def iterate(self):
        if self.sandwich_mode:
            return self.sandwich()
        return self.baguette()

    def get_single_reader(self,number):
        if self.sandwich_mode:
            return self.engine(self.topology,self.trajectory[number],
                              number=number,
                              window=self.window)
        return self.engine(self.topology,self.trajectory,
                          number=0,
                          window=self.window)


class ReaderTraj(object):
    def __init__(self,topology,trajectory,number=None,window=None):
        '''
        :param str topology:  Topology file name.
        :param list trajectory: Trajectory file name.
        :param int number: Number of trajectory file.
        :param Window window: Frames window to read.
        '''

        self.topology = topology
        self.trajectory = trajectory
        if not isinstance(trajectory,list):
            self.trajectory = [trajectory]

        self.number = number
        self.window = window

        self.real_frame_nr = None

        self.trajectory_object = self.open_trajectory()

    def __repr__(self):
        if len(self.trajectory) == 1:
            return "ReaderTraj(%s,%s,%d,%r)" % (self.topology,self.trajectory[0],self.number,self.window)
        return "ReaderTraj(%s,%s,%d,%r)" % (self.topology,"[%s]" % (','.join(self.trajectory)),self.number,self.window)

    def open_trajectory(self):
        # should return any object that can be further used to parse trajectory
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def iterate_over_frames(self):
        for frame,real_frame in enumerate(self.window.range()):
            self.set_real_frame(real_frame)
            yield frame

    def set_real_frame(self,real_frame):
        self.real_frame_nr = real_frame
        # sets real frame
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")


    def parse_selection(self, selection):
        # should return selection
        # selection resolution should be atoms
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
        return [mda.Universe(self.topology_file_name,
                            self.trajectory_file_name,
                            topology_format=topology,
                            format=trajectory)]

    def set_real_frame(self,real_frame):
        self.real_frame_nr = real_frame
        self.trajectory_object.trajectory[real_frame]


s