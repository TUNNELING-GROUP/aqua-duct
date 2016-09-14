from unittest import TestCase
import numpy as np
from traces import LinearizeRecursive


class TestLinearizeRecursive(TestCase):
    #jeszcze nie wiem jak to przetestowac
    def test_here(self):
        test_case=np.array(((2, 3, 4),(2, 3, 6), (2, 3, 10), (2,4,12)))
        res=LinearizeRecursive()

        out=[1,3]
        self.assertEqual(out,res.here(test_case))
