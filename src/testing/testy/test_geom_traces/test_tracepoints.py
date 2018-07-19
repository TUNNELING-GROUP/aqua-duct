from unittest import TestCase

import numpy as np

from aquaduct.geom.traces import tracepoints


class TestTracepoints(TestCase):
    def test_traces1p(self):
        A = np.array([0, 2, 0])
        B = np.array([2, 2, 0])
        self.assertTrue((tracepoints(A, B, 1) == np.array([[1, 2, 0]])).all())

    def test_traces1p_next(self):
        A = np.array([7, 2, 0])
        B = np.array([3, 5, 0])
        self.assertTrue((tracepoints(A, B, 1) == np.array([[5, 3.5, 0]])).all())

    def test_traces1p_more(self):
        A = np.array([0, 2, 0])
        B = np.array([2, 2, 0])
        self.assertTrue((tracepoints(A, B, 3) == np.array([[0.5, 2, 0], [1, 2, 0], [1.5, 2, 0]])).all())
