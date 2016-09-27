from unittest import TestCase
from aqueduct.geom.traces import vectors_angle
from math import pi


class TestVectors_angle(TestCase):
    def test_Vectors_angle_right(self):
        A = (1,1,0)
        B = (0,0,1)
        case = vectors_angle(A,B)
        self.assertEqual(case,pi/2)

    def test_Vectors_angle_acute(self):
        A = (1, 1, 0)
        B = (1, 0, 0)
        case = vectors_angle(A, B)
        self.assertAlmostEqual(case, pi / 4, 6)

    def test_Vectors_angle_0(self):
        A = (0, 0, 1)
        B = (0, 0, 2)
        case = vectors_angle(A, B)
        self.assertEqual(case, 0)

    # def test_Vectors_angle_acute(self):
    #     A = (1, 1)
    #     B = (1, 0)
    #     C = float('nan')
    #     self.assertRaises(ValueError,vectors_angle(A,C))
