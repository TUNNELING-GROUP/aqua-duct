from math import sqrt
from unittest import TestCase

import numpy as np

from aquaduct.geom.traces import diff


class TestDiff(TestCase):
    def test_diff(self):
        A = [0, 2, 0]
        B = [2, 2, 0]
        data = np.array([A, B])
        self.assertEqual(diff(data), 2)

    def test_diff_val(self):
        A = [0, 0, 0]
        B = [2, 2, 0]
        data = np.array([A, B])
        self.assertEqual(diff(data), 2 * sqrt(2))

    def test_diff_4d(self):
        A = [0, 0, 0, 0]
        B = [2, 2, 0, 0]
        data = np.array([A, B])
        self.assertRaises(TypeError, diff(data))
