from unittest import TestCase
from traces import vectors_angle_anorm,vector_norm
import numpy as np

class TestVectors_angle_anorm(TestCase):
    def test_output_value(self):
        A = (0, 2, 2)
        B = (2, 2, 2)
        Anorm = vector_norm(A)
        outp = vectors_angle_anorm(A, B, Anorm)
        self.assertTrue(isinstance(outp,float))

    def test_vectors_angle_anorm(self):
        A=(0,2,2)
        B=(2,2,2)
        Anorm=vector_norm(A)
        outp=vectors_angle_anorm(A,B,Anorm)
        self.assertEqual(outp,np.dot((1,2,2),(3,4,1)))
