from unittest import TestCase
from traces import TrianlgeLinearize


class TestTrianlgeLinearize(TestCase):

    def test_is_linear(self):
        A=(0,2,0)
        B=(2,2,0)
        C=(4,2,0)
        test_case=TrianlgeLinearize(0.05)
        self.assertTrue(test_case.is_linear([A,B,C]))

    def test_false_linear(self):
        A = (0, 2, 0)
        B = (2, 2, 2)
        C = (4, 2, 0)
        test_case = TrianlgeLinearize(0.05)
        self.assertFalse(test_case.is_linear([A, B, C]))

    def test_low_threshold(self):
        A = (0, 2, 0)
        B = (2, 2, 0.05)
        C = (4, 2, 0)
        test_case = TrianlgeLinearize(0.15)
        self.assertTrue(test_case.is_linear([A, B, C]))
    #check the optimal threshold
    #chceck list of list istead of tuple
    #check the numpy.array() input