 
# -*- coding: utf-8 -*-
import unittest
from unittest import TestCase
from aquaduct.geom.pca import Polarize
import numpy as np

class TestPolarize(TestCase):
    def test_polarize_undo(self):
        test_arr = np.random.random((10,3))
        test_cntr = np.random.random(3)
        P = Polarize(center=test_cntr)
        P.build(test_arr)
        test_arr_p = P(test_arr)
        Pp = P.undo(test_arr_p) - test_arr
        [self.assertAlmostEqual(x, 0, 7) for x in Pp.ravel()]
        
    def test_polarize_in(self):
        P = Polarize(center=np.random.random(3))
        self.assertRaises(TypeError, P, 'cupcakes')
        self.assertRaises(TypeError, Polarize, 'cupcakes')
        self.assertRaises(TypeError, Polarize, np.random.random(2))
        self.assertRaises(TypeError, Polarize, np.random.random(4))
        self.assertRaises(ValueError, P, np.random.random((1, 2)))
        self.assertRaises(ValueError, P, np.random.random((1, 4)))
        
    def test_polarize_undo_point(self):
        test_cntr = np.random.random(3)
        test_arr = np.random.random((1,3))
        P = Polarize(center=test_cntr)
        P.build(test_arr)
        test_arr_p = P(test_arr)
        Pp = P.undo(test_arr_p) - test_arr
        [self.assertAlmostEqual(x, 0, 7) for x in Pp.ravel()]



           
if __name__ == '__main__':
    unittest.main()
        
