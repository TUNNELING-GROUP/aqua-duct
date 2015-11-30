'''
Created on Nov 19, 2015

@author: tljm
'''

from aquarium import tests as aqtests

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


if __name__ == "__main__":
    topology = aqtests.get("simple.prmtop")
    trajectory = aqtests.get("simple.nc")
    reader = ReadAmberNetCDF(topology,trajectory)
    print reader.topology_file_name
    print reader.trajectory_file_name
    reader.open_trajectory()
    