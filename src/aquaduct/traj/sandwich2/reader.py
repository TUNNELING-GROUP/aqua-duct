# -*- coding: utf-8 -*-


from abc import ABC, abstractmethod


class BaseReader(ABC):

    def __init__(self,topology,trajectory,
                 frames=None,
                 number=None):
        self.topology = topology
        self.trajectory = trajectory
        self.frames = frames
        self.number = number

        super().__init__()

        self.trajectory_object = self.open_trajectory()



    @abstractmethod
    def open_trajectory(self):
        pass

    @abstractmethod
    def close_trajectory(self):
        pass

    @abstractmethod
    def physical_number_of_frames(self):
        # total number of frames in the supplied trajectory
        pass


    def __del__(self):
        self.close_trajectory()
