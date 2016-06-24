from unittest import TestCase


import numpy as np

from aqueduct.geom.smooth import GeneralWindow


class TestGeneralWindow(TestCase):
    def test_max_window_at_pos(self):
        # TODO: this test shows how simmple is the method - reimplement it!
        for size in range(10,100):
            for pos in range(size/2):
                lo,hi = GeneralWindow.max_window_at_pos(pos,size)
                self.assertEqual(lo,0)
                self.assertEqual(hi,pos*2+1)
            for pos in range(size / 2,size):
                lo, hi = GeneralWindow.max_window_at_pos(pos, size)
                self.assertEqual(hi, size)
                self.assertEqual(lo,pos*2+1-size)


    def test_check_bounds_at_max_window_at_pos(self):
        for size in range(10, 30):
            for lb in range(size):
                for ub in range(lb,size):
                    for pos in range(size / 2):
                        lo,hi = GeneralWindow().check_bounds_at_max_window_at_pos(lb,ub,pos,size)
                        self.assertEqual(lo,0)
                        self.assertEqual(hi, pos+1,'size %d, lb %d, ub %d, pos %d, lo %d, hi %d' % (size,lb,ub,pos,lo,hi))
                    for pos in range(size / 2, size):
                        lo,hi = GeneralWindow().check_bounds_at_max_window_at_pos(lb,ub,pos,size)
                        self.assertEqual(lo,0)

