from unittest import TestCase

from .resources import get

from aquaduct.traj.sandwich2 import MetaReader, available_engines

class TestMetaReader(TestCase):

    def setUp(self):
        self.top = [get('small_top.pdb')]
        self.traj = [get('small_1.dcd'),
                     get('small_2.dcd')]

        self.mrs = []
        for engine in available_engines:
            self.mrs.append(MetaReader(self.top,self.traj,engine=engine))
            self.mrs.append(MetaReader(self.top,self.traj,engine=engine,
                                       mode='sandwich'))

        self.expect_number_of_readers = [1,2] * len(available_engines)

        self.expect_physical_number_of_frames = [4000,4000] * len(available_engines)


    def test_number_of_readers(self):
        for mr,e in zip(self.mrs,self.expect_number_of_readers):
            self.assertEqual(mr.number_of_readers(),e)

    def test_physical_number_of_frames(self):
        for mr,e in zip(self.mrs,self.expect_physical_number_of_frames):
            self.assertEqual(mr.physical_number_of_frames(),e)
        pass

