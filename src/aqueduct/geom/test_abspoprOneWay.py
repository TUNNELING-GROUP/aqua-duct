from unittest import TestCase
from traces import abspoprOneWay
import numpy as np
import types
from itertools import islice

class TestAbspoprOneWay(TestCase):
    def test_here0(self):
        A = [0, 2, 0]
        B = [2, 2, 0]
        C = [4, 2, 0]
        data = np.array([A, B, C])
        test_case = abspoprOneWay(0.5)
        self.assertTrue(isinstance(test_case.here(data),types.GeneratorType))

    def test_here1(self):
        A = [0, 2, 0]
        B = [2, 2, 0]
        C = [4, 2, 0]
        out = (0,2)
        data = np.array([A, B, C])
        test_case = abspoprOneWay(0.5)
        self.assertEqual(tuple(test_case.here(data)),out)

    def test_here2(self):
        A = [0, 2, 0]
        B = [2, 2, 0]
        C = [4, 7, 5]
        D = [5, 7, 1]
        E = [6, 7, 2]
        F = [7, 7, 2]
        G = [8, 7, 2]
        H = [9, 7, 2]
        I = [9, 8, 3]
        out = (0,1,2,3,4,7,8) #4,8
        data = np.array([A, B, C, D, E, F, G, H, I])
        test_case = abspoprOneWay(0)
        self.assertEqual(tuple(test_case.here(data)), out)

    def test_here3(self):
        A = [0, 2, 0]
        B = [2, 2, 0]
        C = [4, 2, 0]
        D = [6, 2, 0]
        E = [4, 7, 0]
        F = [2, 3, 2]
        G = [4, 3, 1]
        test_case = abspoprOneWay(0.1)
        out = (0,3,4,5,6)
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
        out = (0, 5,6,7,8)
        data = np.array([A, B, C, D, E, F, G, H, I])
        test_case = abspoprOneWay(0.00)
        self.assertEqual(tuple(test_case.here(data)), out)


    def test_here5(self):
        A = [0, 2, 0]
        B = [2, 5, 0]
        C = [4, 2, 0]
        D = [5, 5, 0]
        E = [6, 2, 0]
        F = [7, 2, 0]
        G = [8, 2, 0]
        H = [9, 2, 0]
        I = [9, 8, 3]
        J = [7, 5, 6]
        K = [8, 2, 9]
        L = [6, 4, 4]
        M = [7, 4, 4]
        N = [8, 4, 4]
        O = [9, 4, 4]
        out = (0,1,2,3,4,7,8,9,10,11,14)
        data = np.array([A, B, C, D, E, F, G, H, I,J,K,L,M,N,O])
        test_case = abspoprOneWay(0.5)
        self.assertEqual(tuple(test_case.here(data)), out)