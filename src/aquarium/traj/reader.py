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

        self.trajectory_object = self.open_trajectory()
        self.set_current_frame(0) # ensure that by default it starts from frame 0

    def open_trajectory(self):
        # should return any object that can be further used to parse trajectory
        raise NotImplementedError()
    
    @property
    def number_of_frames(self):
        # should return number of frames
        raise NotImplementedError()

    def set_current_frame(self,frame):
        # should set current frame
        raise NotImplementedError()

    def next_frame(self):
        # should set next frame or raise StopIteration
        raise NotImplementedError()

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

    def parse_selection(self,selection):
        # should parse and return selection object
        # is this object updated automatically depends on the particular implementation of Reader class
        raise NotImplementedError()

class ReadAmberNetCDFviaMDA(Reader):

    def open_trajectory(self):
        return mda.Universe(self.topology_file_name,
                            self.trajectory_file_name,
                            topology_format="prmtop",
                            format="nc")

    @property
    def number_of_frames(self):
        # should return number of frames
        return len(self.trajectory_object.trajectory)

    def set_current_frame(self,frame):
        self.trajectory_object.trajectory[frame]

    def next_frame(self):
        try:
            self.trajectory_object.trajectory.next()
        except:
            raise StopIteration

    def parse_selection(self,selection):
        return self.trajectory_object.select_atoms(selection)

