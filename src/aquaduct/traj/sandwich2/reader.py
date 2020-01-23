# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2018-2019  Tomasz Magdziarz, Michał Banas <info@aquaduct.pl>
# Copyright (C) 2020  Tomasz Magdziarz <info@aquaduct.pl>
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


from aquaduct.traj.sandwich2.raw import ENGINE_MDA, open_raw

################################################################################

class ReaderTraj(object):

    def __init__(self, topology=None, trajectory=[],
                 window=slice(None)):
        self.topology = topology
        self.trajectory = trajectory
        self.window = window

        self.trajectory_object = self.open()

    def open(self):
        # should return any object that can be further used to parse trajectory
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def close(self):
        # should close trajectory reader in self.trajectory_object
        # WARNING: This method has to be carefully implemented because it is used by __del__ and
        #          should not emmit any error messages. This is subjet of change.
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def __del__(self):
        return self.close()


################################################################################

class ReaderTrajViaMDA(ReaderTraj):

    def open(self):
        return open_raw(self.topology,self.trajectory, ENGINE_MDA)

    def close(self):
        if hasattr(self, 'trajectory_object'):
            if hasattr(self.trajectory_object, 'trajectory'):
                if hasattr(self.trajectory_object.trajectory, 'close'):
                    self.trajectory_object.trajectory.close()


################################################################################
