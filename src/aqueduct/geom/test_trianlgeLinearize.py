from unittest import TestCase
from traces import TriangleLinearize
import numpy as np
import math

class TestTrianlgeLinearize(TestCase):

    def test_is_linear(self):
        A=(0,2,0)
        B=(2,2,0)
        C=(4,2,0)
        test_case=TriangleLinearize(0.05)
        self.assertTrue(test_case.is_linear([A,B,C]))

    def test_false_linear(self):
        A = (0, 2, 0)
        B = (2, 2, 2)
        C = (4, 2, 0)
        test_case = TriangleLinearize(0.05)
        self.assertFalse(test_case.is_linear([A, B, C]))

    def test_low_threshold(self):
        A = (0, 2, 0)
        B = (2, 2, 0.05)
        C = (4, 2, 0)
        test_case = TriangleLinearize(0.15)
        self.assertTrue(test_case.is_linear([A, B, C]))

    def test_list(self):
        A = [0, 2, 0]
        B = [2, 2, 0.1]
        C = [4, 2, 0]
        test_case= TriangleLinearize(0.5)
        self.assertTrue(test_case.is_linear([A, B, C]))

    def test_numpyarr(self):
        A = [0, 2, 0]
        B = [2, 2, 0.1]
        C = [4, 2, 0]
        A, B, C = np.array([A, B, C])
        test_case = TriangleLinearize(0.5)
        self.assertTrue(test_case.is_linear([A, B, C]))

    def test_sillyvalue(self):
        #dla dwoch punktow tez dziala
        A = [0, 2, 0]
        B = [2, 2, 0.1]
        C = float('nan')
        test_case = TriangleLinearize(0.5)
        #self.assertTrue(math.isnan(C))
        self.assertRaises(TypeError,test_case.is_linear([A, B, C]))