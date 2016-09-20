from unittest import TestCase
from traces import VectorLinearize
import numpy as np

class TestVectorLinearize(TestCase):

    def test_threshold_type(self):
        a=0
        obj=VectorLinearize(a)
        self.assertTrue(isinstance(obj,VectorLinearize))

    # todo: threshold moze byc ujemny?
    # def test_threshold_negative(self):
    #     a = -100
    #     self.assertRaises(ValueError, VectorLinearize(-10))

    def test_threshold_abovezero(self):
        a = 0
        A = [0, 2, 0]
        B = [2, 2, 0]
        C = [4, 2, 0]
        D = [5, 2, 0]
        E = [6, 2, 0]
        data=[A,B,C,D,E]
        typed=np.array(data)
        obj = VectorLinearize(a)
        self.assertTrue(obj.is_linear(typed)==True)

    def test_threshold_CORE(self):
        a = 0
        A = [0, 2, 0]
        B = [2, 2, 0]
        C = [4, 2, 0]
        D = [5, 2, 0]
        E = [6, 2, 0]
        data = [A, B, C, D,E ]
        typed = np.array(data)
        obj = VectorLinearize(a)
        self.assertTrue(obj.is_linear_core(typed) == True)

    def test_threshold_false(self):
        a = 0
        A = [0, 2, 0]
        B = [2, 2, 0]
        C = [4, 8, 0]
        D = [5, 2, 0]
        E = [6, 2, 0]
        data=[A,B,C,D,E]
        typed=np.array(data)
        obj = VectorLinearize(a)
        self.assertTrue(obj.is_linear(typed)==False)

    def test_threshold_false_core(self):
        a = 0
        A = [0, 2, 0]
        B = [2, 2, 0]
        C = [4, 8, 0]
        D = [5, 2, 0]
        E = [6, 2, 0]
        data = [A, B, C, D,E ]
        typed = np.array(data)
        obj = VectorLinearize(a)
        self.assertTrue(obj.is_linear_core(typed) == False)

    #TODO: czy powinien zwracac true dla 2 wartosci?
    def test_2val(self):
        a=0
        A = [0, 2, 0]
        B = [2, 2, 0]
        data = [A, B]
        typed = np.array(data)
        obj = VectorLinearize(a)
        self.assertFalse(obj.is_linear(typed)==True)
g
