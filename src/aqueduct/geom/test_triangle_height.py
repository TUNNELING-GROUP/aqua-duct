# -*- coding: utf-8 -*-
from unittest import TestCase
from traces import triangle_height
import math

class TestTriangle_height(TestCase):
    def test_OutType(self):
        A = 1
        B = 2
        C = 3
        case = triangle_height(A, B, C)
        self.assertTrue(isinstance(case, float))
    def test_2dim(self):
        A=0,2
        B=2,2
        C=2,4
        case=triangle_height(A,B,C)
        self.assertEquals(case, 2)
    def test_3dim(self):
        A = 0, 2, 0
        B = 2, 2, 2
        C = 2, 4, 2
        case = triangle_height(A, B, C)
        self.assertEquals(case,2.449489743)