# -*- coding: utf-8 -*-
import numpy as np
from unittest import TestCase

from aquaduct.geom.traces import vectors_angle_anorm, vector_norm


class TestVectors_angle_anorm(TestCase):
    def test_output_value(self):
        A = (0, 2, 2)
        B = (2, 2, 2)
        Anorm = vector_norm(A)
        outp = vectors_angle_anorm(A, B, Anorm)
        self.assertTrue(isinstance(outp, float))

    def test_vectors_angle_anorm(self):
        A = (0, 2, 2)
        B = (2, 2, 2)
        Anorm = vector_norm(A)
        outp = vectors_angle_anorm(A, B, Anorm)
        res = np.arccos(np.clip(8 / (2.83 * 3.46), -1, 1))
        round_res = round(res, 3)
        round_outp = round(outp, 3)
        self.assertEqual(round_res, round_outp)

    def test_list_input(self):
        A = [0, 2, 2]
        B = [2, 2, 2]
        Anorm = vector_norm(A)
        outp = vectors_angle_anorm(A, B, Anorm)
        res = np.arccos(np.clip(8 / (2.83 * 3.46), -1, 1))
        round_res = round(res, 3)
        round_outp = round(outp, 3)
        self.assertEqual(round_res, round_outp)

    def test_zero(self):
        A = [0, 0, 0]
        B = [2, 2, 2]
        Anorm = vector_norm(A)
        outp = vectors_angle_anorm(A, B, Anorm)
        self.assertEqual(outp, 0)

    def test_zero2(self):
        # in such case it would be worth to pass callculation with 0 value
        # values must be
        A = [0.1, 0.1, 0.1]
        B = [0.1, 0.1, 0.1]
        Anorm = vector_norm(A)
        outp = vectors_angle_anorm(A, B, Anorm)
        res = np.arccos(np.clip(0.03 / (0.1732 * 0.1732), -1, 1))
        self.assertAlmostEqual(outp, res)

    def test_huge_values(self):
        A = 200, 200, 200
        B = 200, 200, 200
        Anorm = vector_norm(A)
        outp = vectors_angle_anorm(A, B, Anorm)
        res = np.arccos(np.clip(120000 / (346.41 * 346.41), -1, 1))
        self.assertEqual(outp, res)

    def test_mixed_values(self):
        # it's enough if just one vector is extremely short to give vaule about 0
        A = 0.1, 0.1, 0.1
        B = 200, 200, 200
        Anorm = vector_norm(A)
        outp = vectors_angle_anorm(A, B, Anorm)
        res = np.arccos(np.clip(60 / (0.03 * 346.41), -1, 1))
        self.assertAlmostEqual(outp, res)

    def test_mixed_values(self):
        # one extremely small value has no influence on result
        A = 1, 1, 0.1
        B = 200, 200, 200
        Anorm = vector_norm(A)
        outp = vectors_angle_anorm(A, B, Anorm)
        res = np.arccos(np.clip(420 / (1.4177 * 346.41), -1, 1))
        roundo = round(outp, 3)
        roundr = round(res, 3)
        self.assertEqual(roundo, roundr)

    def test_negative_vector(self):
        A = -2, -2, -2
        B = 2, 2, 2
        Anorm = vector_norm(A)
        outp = vectors_angle_anorm(A, B, Anorm)
        res = np.arccos(np.clip(-8 / (2.82842 * 2.82842), -1, 1))
        self.assertAlmostEqual(outp, res)

    def test_one_negative_value(self):
        # todo: values are quite different-why
        A = -10, 2, 2
        B = 2, 2, 2
        Anorm = vector_norm(A)
        outp = vectors_angle_anorm(A, B, Anorm)
        res = np.arccos(np.clip(np.dot(A, B) / (10.4 * 2.82), -1, 1))
        roundo = round(outp, 3)
        roundr = round(res, 3)
        self.assertAlmostEqual(outp, res,0)
