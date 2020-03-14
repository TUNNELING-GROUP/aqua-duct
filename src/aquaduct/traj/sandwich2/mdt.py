# -*- coding: utf-8 -*-

from aquaduct import logger


from aquaduct.traj.sandwich2.reader import BaseReader

from os.path import splitext
import re



import mdtraj as md



################################################################################
# raw

available_formats = {re.compile('(dcd|DCD)'): 'LAMMPS',
                     re.compile('(pdb|PDB)'): 'pdb',
                     re.compile('(xtc|XTC)'): 'XTC'}

def open_raw(topology, trajectory):
    topology_ext = splitext(topology)[1][1:]
    for afk in list(available_formats.keys()):
        if afk.match(topology_ext):
            topology_ext = available_formats[afk]
            break
    trajectory_ext = splitext(trajectory[0])[1][1:]
    for afk in list(available_formats.keys()):
        if afk.match(trajectory_ext):
            trajectory_ext = available_formats[afk]
            break

    return md.load(trajectory,top=topology)

################################################################################


class Reader(BaseReader):

    def open_trajectory(self):
        # returns raw trajectory objet to be interpreted by this class
        return open_raw(self.topology,self.trajectory)


    def close_trajectory(self):
        if hasattr(self, 'trajectory_object'):
            if hasattr(self.trajectory_object, 'trajectory'):
                if hasattr(self.trajectory_object.trajectory, 'close'):
                    self.trajectory_object.trajectory.close()

    def physical_number_of_frames(self):
        return len(self.trajectory_object.trajectory)
