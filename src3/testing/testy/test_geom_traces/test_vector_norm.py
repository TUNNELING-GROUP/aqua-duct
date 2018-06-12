# -*- coding: utf-8 -*-
from unittest import TestCase
from aquaduct.geom.traces import vector_norm
import math

class TestVector_norm(TestCase):
    def test_vector_norm(self):
        test_val=(5)
        case1=vector_norm(test_val)
        self.assertEqual(case1,5)

    def test_tuple_norm(self):
        #input can be a tuple...
        test_val=(4,3)
        case=vector_norm(test_val)
        self.assertEqual(case,5)

    def test_list(self):
        #...or it can be a list
        test_val=[4,3]
        case=vector_norm(test_val)
        self.assertEqual(case,5)
    def test_large_triangle(self):
        test_val= (9,12)
        case=vector_norm(test_val)
        self.assertEqual(case,15)

    def test_3d(self):
        #calculate 3 dimensional vector length too
        test_val=(2,3,4)
        case=vector_norm(test_val)
        self.assertAlmostEqual(case, 5.3851648)

    def test_type(self):
        test_val = (5)
        case1 = vector_norm(test_val)
        self.assertTrue(case1, float)

    def test_tupletuple(self):
        test_val = ((0,0,0),(0,2,2))
        case1 = vector_norm(test_val)
        self.assertEqual(case1,2*math.sqrt(2) )
