from unittest import TestCase
from aquaduct.geom.traces import length_step_std
import numpy as np
import math


class TestLength_step_std(TestCase):
    def test_length_step_std(self):
        A = [0, 0, 0]
        B = [2, 0, 0]
        data = np.array([A, B])
        d = length_step_std(data)
        self.assertEqual((2, 2, 0), d)

    def test_length_step_std1(self):
        A = [0, 0, 0]
        B = [2, 0, 0]
        C = [4, 0, 0]
        D = [6, 0, 0]
        trace = np.array([A, B, C, D])
        d = length_step_std(trace)
        self.assertEqual((6, 2, 0), d)

    def test_length_step_std2(self):
        A = [0, 0, 0]
        B = [2, 2, 2]
        C = [4, 4, 4]
        D = [6, 6, 6]
        trace = np.array([A, B, C, D])
        d = length_step_std(trace)
        dt = tuple([round(point) for point in d])
        self.assertEqual((round(3 * 2 * math.sqrt(3)), round(2 * math.sqrt(3)), 0), dt)
