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
        out = (0,1,2)
        data = np.array([A, B, C])
        test_case = absOneWay(0.5)
        self.assertEqual(tuple(islice(test_case.here(data),1)),out)
