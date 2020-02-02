# -*- coding: utf-8 -*-



from itertools import zip_longest



ENGINE_MDA = 'mda'
available_engines = [ENGINE_MDA]



################################################################################
# Trajectory data

class MetaReader(object):

    baguette = 'baguette'
    sandwich = 'sandwich'
    shortbread = 'shortbread'

    def __init__(self, topology=[], trajectory=[],
                 mode='baguette',
                 window=slice(None),
                 threads=1,
                 engine=ENGINE_MDA):

        self.topology = topology
        self.trajectory = trajectory
        self.mode = mode
        self.window = window
        self.threads = threads
        self.engine = engine
        # engine initialization pending
        # would be nice to know number of frames

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
    # this is not a real reader, it can be safely passed becuse it is a very simple object

    def __init__(self, topology, trajectory,
                 frames=None,
                 number=None,
                 engine=None):
        self.topology = topology
        self.trajectory = trajectory
        self.frames = frames
        self.number = number
        self.engine = engine

    def open(self):
        if self.engine == ENGINE_MDA:
            from aquaduct.traj.sandwich2.mda import Reader
        # return reader for current engine
        return Reader(self.topology,self.trajectory,
                      frames=self.frames,
                      number=self.number)


