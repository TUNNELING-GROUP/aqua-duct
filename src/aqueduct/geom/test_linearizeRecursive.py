from unittest import TestCase
import numpy as np
from traces import LinearizeRecursive
from traces import TriangleLinearize


class TestLinearizeRecursive(TestCase):
    #jeszcze nie wiem jak to przetestowac
    def test_here(self):
        A = [0, 2, 0]
        B = [2, 2, 0]
        C = [4, 2, 0]
        test_case =  LinearizeRecursive(0.5)
        out=[0,1,2]
        data = np.array([A,B, C])
        self.assertEqual(out,test_case.here(data))
# todo: "AttributeError: 'LinearizeRecursive' object has no attribute 'is_linear'
    def test_here_nonlinear(self):
        A = [0, 0, 0]
        B = [2, 0, 0]
        C = [4, 0, 0]
        D = [8, 6, 8]
        E = [16, 5, 0]
        test_case =  LinearizeRecursive(0.5)
        out=[0,1,2]
        data = np.array([A, B, C, D, E])
        self.assertEqual(out,test_case.here(data))
