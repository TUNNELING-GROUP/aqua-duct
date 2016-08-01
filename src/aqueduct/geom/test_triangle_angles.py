from unittest import TestCase
from traces import triangle_angles
import math
class TestTriangle_angles(TestCase):
    def test_type(self):
        A = 1
        B = 2
        C = 3
        case = triangle_angles(A, B, C)
        self.assertTrue(isinstance(case,list))
    def test_triangle_angles_1d(self):
        #1 dimensional case
        A=1
        B=2
        C=3
        case = triangle_angles(A,B,C)
        self.assertEqual(case,[0.0,0.0,math.pi])

    def test_triangle_angles_2d(self):
        #two dimensional
        A=(0,2)
        B=(2,2)
        C=(2,4)
        case=triangle_angles(A,B,C)
        #output is not rounded assertEqual to raw output returns false
        myRoundedList = [round(elem, 2) for elem in case]
        result=[math.pi/4,math.pi/4,math.pi/2]
        myRoundedResult=[round(elem,2) for elem in result]
        self.assertEqual(myRoundedList,myRoundedResult)

    def test_triangle_angles2_d(self):
        # two dimensional
        A = (0, 2)
        B = (2, 2)
        C = (2, 3)
        case = triangle_angles(A, B, C)
        myRoundedList = [round(elem,1) for elem in case]
        result = [0.15*math.pi, 0.35*math.pi, math.pi / 2]
        myRoundedResult = [round(elem, 1) for elem in result]
        self.assertEqual(myRoundedList, myRoundedResult)

    def test_triangle_angles3d(self):
        A = (0, 2, 0)
        B = (2, 2, 0)
        C = (2, 3, 0)
        case = triangle_angles(A, B, C)
        myRoundedList = [round(elem, 1) for elem in case]
        result = [0.15 * math.pi, 0.35 * math.pi, math.pi / 2]
        myRoundedResult = [round(elem, 1) for elem in result]
        self.assertEqual(myRoundedList, myRoundedResult)

    def test_triangle_angles_list(self):
        A = [0, 2, 0]
        B = [2, 2, 0]
        C = [2, 3, 0]
        case = triangle_angles(A, B, C)
        myRoundedList = [round(elem, 1) for elem in case]
        result = [0.15 * math.pi, 0.35 * math.pi, math.pi / 2]
        myRoundedResult = [round(elem, 1) for elem in result]
        self.assertEqual(myRoundedList, myRoundedResult)