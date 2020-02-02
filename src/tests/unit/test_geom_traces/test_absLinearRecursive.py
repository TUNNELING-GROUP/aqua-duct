from unittest import TestCase
from aquaduct.geom.traces import LinearizeRecursiveTriangle
import numpy as np


class TestAbsLinearRecursive(TestCase):
    def test_here(self):
        A = [0, 2, 0]
        B = [2, 2, 0]
        C = [4, 2, 0]
        test_case = LinearizeRecursiveTriangle(0.5)
        out = [0, 1, 2]
        data = np.array([A, B, C])
        self.assertEqual(out, test_case.here(data))

    def test_here2(self):
        # DIDN'T PASS
        A = [0, 2, 0]
        B = [2, 2, 0]
        C = [4, 2, 0]
        D = [3, 6, 0]
        E = [4, 7, 0]
        F = [2, 3, 2]
        G = [4, 3, 1]
        test_case = LinearizeRecursiveTriangle(0.5)
        out = [0, 2, 3, 4, 5, 6]
        data = np.array([A, B, C, D, E, F, G])
        self.assertEqual(out, test_case.here(data))

    def test_here3(self):
        # DIDN'T PASS
        A = [0, 2, 0]
        B = [2, 2, 0]
        C = [4, 2, 0]
        C1 = [6, 2, 0]
        D = [3, 6, 0]
        E = [4, 7, 0]
        F = [2, 3, 2]
        G = [4, 3, 1]
        test_case = LinearizeRecursiveTriangle(0.05)
        out = [0, 3, 4, 5, 6, 7]
        data = np.array([A, B, C, C1, D, E, F, G])
        self.assertEqual(out, test_case.here(data))

    def test_here4(self):
        A = [0, 2, 0]
        B = [2, 2, 0]
        B1 = [3, 2, 0]
        C = [4, 2, 0]
        D = [3, 6, 0]
        test_case = LinearizeRecursiveTriangle(0.5)
        out = [0, 2, 3, 4]
        data = np.array([A, B, B1, C, D])
        self.assertEqual(out, test_case.here(data))

    def test_here5(self):
        A = [0, 2, 0]
        B = [2, 2, 0]
        B1 = [3, 2, 0]
        C = [4, 2, 0]
        D = [3, 6, 0]
        D1 = [3, 8, 0]
        D2 = [3, 10, 0]
        test_case = LinearizeRecursiveTriangle(0.5)
        out = [0, 2, 3, 4, 5, 6]
        data = np.array([A, B, B1, C, D, D1, D2])
        self.assertEqual(out, test_case.here(data))
