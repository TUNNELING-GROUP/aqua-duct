'''
Created on Nov 19, 2015

@author: tljm
'''

import MDAnalysis as mda

class ReadAmberNC(object):
    def __init__(self,topology,trajectory):
        assert isinstance(topology, string)
        assert isinstance(trajectory, string)
        
        self.topology_file_name = topology
        self.trajectory_file_name = trajectory