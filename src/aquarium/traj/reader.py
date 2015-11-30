'''
Created on Nov 19, 2015

@author: tljm
'''

import MDAnalysis as mda

class Reader(object):
    
    def __init__(self,topology,trajectory):
        assert isinstance(topology, str)
        assert isinstance(trajectory, str)
        
        self.topology_file_name = topology
        self.trajectory_file_name = trajectory
    
    def open_trajectory(self):
        raise NotImplementedError()
    
        
class ReadAmberNetCDF(Reader):
    pass




reader = ReadAmberNetCDF("topology","trajectory")
reader.open_trajectory()