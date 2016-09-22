from unittest import TestCase
from traces import vectors_angle
from math import pi


class TestVectors_angle(TestCase):
    def test_Vectors_angle(self):
        A = (1,1,0)
        B = (0,0,1)
        case = vectors_angle(A,B)
        self.assertEqual(case,pi/2)
    def test_Vectors_angle_3val(self):
        A = (1,1,0)
        B = (0,0,1)
        C = (2,1,0)
        self.assertRaises()
