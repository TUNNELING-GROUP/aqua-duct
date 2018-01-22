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

class Window(object):
    def __init__(self,start,stop,step):
        self.start = start
        self.stop = stop
        self.step = step

    def __repr__(self):
        return "Window(%r:%r:%r)" % (self.start,self.stop,self.step)

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
        '''

        self.topology = topology
        self.trajectory = trajectory

        self.window = window
        self.sandwich_mode = sandwich

    def __repr__(self):
        sandwich = ''
        if self.sandwich_mode:
            sandwich = ',sandwich'
        return "Reader(%s,%s,%r%s)" % (self.topology,"[%s]" % (','.join(self.trajectory)),self.window,sandwich)

    def sandwich(self):
        for nr,trajectory in enumerate(self.trajectory):
            yield ReaderTraj(self.topology,trajectory,
                              number=nr,
                              window=self.window)

    def baguette(self):
        yield ReaderTraj(self.topology,self.trajectory,
                          number=0,
                          window=self.window)

    def iterate(self):
        if self.sandwich_mode:
            return self.sandwich()
        return self.baguette()

    def get_single_reader(self,number):
        if self.sandwich_mode:
            return ReaderTraj(self.topology,self.trajectory[number],
                              number=number,
                              window=self.window)
        return ReaderTraj(self.topology,self.trajectory,
                          number=0,
                          window=self.window)


class ReaderTraj(object):
    def __init__(self,topology,trajectory,number=None,window=None):
        '''
        :param str topology:  Topology file name.
        :param trajectory: Trajectory file name.
        :param int number: Number of trajectory file.
        :param Window window: Frames window to read.
        '''

        self.topology = topology
        self.trajectory = trajectory
        if not isinstance(trajectory,list):
            self.trajectory = [trajectory]
        self.number = number

        self.window = window

    def __repr__(self):
        if len(self.trajectory) == 1:
            return "ReaderTraj(%s,%s,%d,%r)" % (self.topology,self.trajectory[0],self.number,self.window)
        return "ReaderTraj(%s,%s,%d,%r)" % (self.topology,"[%s]" % (','.join(self.trajectory)),self.number,self.window)

