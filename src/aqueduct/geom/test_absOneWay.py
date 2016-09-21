from unittest import TestCase
from traces import absOneWay
import numpy as np
import types
from itertools import islice

class TestAbsOneWay(TestCase):
    def test_here0(self):
        A = [0, 2, 0]
        B = [2, 2, 0]
        C = [4, 2, 0]
        data = np.array([A, B, C])
        test_case = absOneWay(0.5)
        self.assertTrue(isinstance(test_case.here(data),types.GeneratorType))

    def test_here1(self):
        A = [0, 2, 0]
        B = [2, 2, 0]
        C = [4, 7, 0]
        out = (0,)
        data = np.array([A, B, C])
        test_case = absOneWay(0.5)
        self.assertEqual(tuple(islice(test_case.here(data),1)),out)

    def test_here2(self):
        A = [0, 2, 0]
        B = [2, 2, 0]
        C = [4, 7, 0]
        D = [5, 7, 1]
        E = [6, 7, 2]
        F = [7, 7, 2]
        G = [8, 7, 2]
        H = [9, 7, 2]
        I = [9, 8, 3]
        out = (0,3,7) #4,8
        data = np.array([A, B, C, D, E, F, G, H, I])
        test_case = absOneWay(0.5)
        self.assertEqual(tuple(test_case.here(data)), out)

    def test_here3(self):
        A = [0, 2, 0]
        B = [2, 2, 0]
        C = [4, 2, 0]
        D = [3, 6, 0]
        E = [4, 7, 0]
        F = [2, 3, 2]
        G = [4, 3, 1]
        test_case = absOneWay(0.1)
        out = (0, 4)
        data = np.array([A, B, C, D, E, F, G])
        self.assertEqual(out, tuple(test_case.here(data)))


    def test_here4(self):
        A = [0, 2, 0]
        B = [2, 2, 0]
        C = [4, 2, 0]
        D = [5, 2, 0]
        E = [6, 2, 0]
        F = [7, 2, 0]
        G = [8, 6, 3]
        H = [9, 7, 2]
        I = [9, 8, 3]
        out = (0, 4)
        data = np.array([A, B, C, D, E, F, G, H, I])
        test_case = absOneWay(0.00)
        self.assertEqual(tuple(test_case.here(data)), out)