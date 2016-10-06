from unittest import TestCase
from aqueduct.geom.traces import tracepoints
import numpy as np


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
        self.assertTrue((tracepoints(A, B, 3) == np.array([[0.5, 2, 0],[1,2,0],[1.5,2,0]])).all())
