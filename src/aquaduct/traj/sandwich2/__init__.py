# -*- coding: utf-8 -*-

"""
Provides frontend tools for accessing trajectory data.
This should be the only point for accessing this kind of data.
"""

from itertools import zip_longest

################################################################################
# List of engines

ENGINE_MDA = 'mda'
"""Code name for MDAnalysis engine."""
ENGINE_MDT = 'mdt'
"""Code name for MDTraj engine."""
available_engines = [ENGINE_MDA, ENGINE_MDT]
"""List of available MD engines."""


################################################################################
# Trajectory data

class MetaReader(object):
    """
    Main reader object. Provides readers according to settings.
    """

    baguette = 'baguette'
    """Code name for standard baguette mode."""
    sandwich = 'sandwich'
    """Code name for sandwich mode."""
    shortbread = 'shortbread'
    """Code name for waterfall mode."""

    def __init__(self, topology=[], trajectory=[],
                 mode='baguette',
                 window=slice(None),
                 threads=1,
                 engine=ENGINE_MDA):
        """
        :param topology: List of topology file names.
        :param trajectory: List of trajectory file names.
        :param frames: Range of frames reader should return. If None all frames are returned.
        :param number: Number of reader.
        :param threads: Number of threads for processing trajectory data.
        :param engine: Name of MD engine.
        """

        self.topology = topology
        self.trajectory = trajectory
        self.mode = mode
        self.window = window
        self.threads = threads
        self.engine = engine

    def number_of_readers(self):
        if self.mode == self.sandwich:
            return len(self.trajectory)
        return self.threads

    def physical_number_of_frames(self):
        nof = 0
        for top,traj in self.iter_top_traj_pairs():
            reader = self.get_reader(top,traj)
            nof += reader.open().physical_number_of_frames()
        return nof

    def get_reader(self,*args,**kwargs):
        return ProtoReader(engine=self.engine, *args, **kwargs)

    def iter_top_traj_pairs(self):
        # iterates over possible topology and trajectory pairs
        if self.mode in [self.sandwich, self.shortbread]:
            for top,traj in zip_longest(self.topology,
                                        self.trajectory,
                                        fillvalue=self.topology[0]):
                yield top,[traj]
        elif self.mode == self.baguette:
            yield self.topology[0],self.trajectory


################################################################################

class ProtoReader(object):
    """
    Provides prototype of Reader object. Its only function is to open a real Reader
    object with :py:meth:`open` method.
    """

    def __init__(self, topology, trajectory,
                 frames=None,
                 number=None,
                 engine=None):
        """
        :param topology: Topology file name.
        :param trajectory: List of trajectory file names.
        :param frames: Range of frames reader should return. If None all frames are returned.
        :param number: Number of reader.
        :param engine: Name of MD engine.
        """
        self.topology = topology
        self.trajectory = trajectory
        self.frames = frames
        self.number = number
        self.engine = engine

    def open(self):
        """
        Checks the engine settings and imports appropriate Reader object.
        This object iitialized with topology, trajectory, frames, and number data is returned.

        :return: Reader object.
        """
        if self.engine == ENGINE_MDA:
            from aquaduct.traj.sandwich2.mda import Reader
        # return reader for current engine
        return Reader(self.topology,self.trajectory,
                      frames=self.frames,
                      number=self.number)


