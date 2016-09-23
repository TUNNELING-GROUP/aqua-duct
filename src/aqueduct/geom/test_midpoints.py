from unittest import TestCase
from traces import midpoints
import numpy as np


class TestMidpoints(TestCase):
    def test_midpoints0(self):
        A = np.array([[1, 1, 1], [3, 3, 3], [5, 5, 5]])
        B = np.array([[2, 2, 2], [3, 3, 3], [5, 5, 5]])
        l1 = (A, B)
        a = tuple(midpoints(l1))
        out=np.array([[1, 1, 1], [3, 3, 3], [5, 5, 5],[3.5, 3.5, 3.5]])
        b=np.array(a[0])
        self.assertTrue(out.all() == b.all())

    def test_midpoints01(self):
        A = np.array([[1, 1, 1], [3, 3, 3], [5, 5, 5]])
        B = np.array([[2, 2, 2], [3, 3, 3], [5, 5, 5]])
        l1 = (A, B)
        a = tuple(midpoints(l1))
        out = np.array([[3.5, 3.5, 3.5],[2, 2, 2], [3, 3, 3], [5, 5, 5]])
        b = np.array(a[1])
        self.assertTrue(out.all() == b.all())

    def test_midpoints1_1(self):
        A = np.array([[0, 0, 1], [0, 0, 2], [0, 0, 3]])
        B = np.array([[0, 0, 2], [0, 1, 3], [1, 0, 5]])
        l1 = (A, B)
        a = tuple(midpoints(l1))
        out=np.array([[0, 0, 1], [0, 0, 2], [0, 0, 3],[0, 0, 2.5]])
        b = np.array(a[0])
        self.assertTrue(out.all() == b.all())

    def test_midpoints1_2(self):
        A = np.array([[0, 0, 1], [0, 0, 2], [0, 0, 3]])
        B = np.array([[0, 0, 2], [0, 1, 3], [1, 0, 5]])
        l1 = (A, B)
        a = tuple(midpoints(l1))
        out = np.array([[0, 0, 2.5],[0, 0, 2], [0, 1, 3], [1, 0, 5]])
        b = np.array(a[1])
        self.assertTrue(out.all() == b.all())

    def test_midpoints2_1(self):
        A = np.array([[0, 0, 1], [0, 0, 2], [3, 5, 3]])
        B = np.array([[5, 7, 5], [0, 1, 3], [1, 0, 5]])
        l1 = (A, B)
        a = tuple(midpoints(l1))
        out=np.array([[0, 0, 1], [0, 0, 2], [3, 5, 3],[4,6,4]])
        b = np.array(a[0])
        self.assertTrue(out.all() == b.all())

    def test_midpoints2_2(self):
        A = np.array([[0, 0, 1], [0, 0, 2], [3, 5, 3]])
        B = np.array([[5, 7, 5], [0, 1, 3], [1, 0, 5]])
        l1 = (A, B)
        a = tuple(midpoints(l1))
        out = np.array([[4,6,4],[5, 7, 5], [0, 1, 3], [1, 0, 5]])
        b = np.array(a[1])
        self.assertTrue(out.all() == b.all())

    def test_midpoints3_1(self):
        A = np.array([[1, 1, 1], [3, 3, 3], [5, 5, 5]])
        B = np.array([[2, 2, 2], [3, 3, 3], [5, 5, 5]])
        C = np.array([[3, 3, 3], [5, 5, 5], [7, 7, 7]])
        l1 = (A, B, C)
        a = tuple(midpoints(l1))
        out = np.array([[1, 1, 1], [3, 3, 3], [5, 5, 5], [3.5, 3.5, 3.5]])
        b = np.array(a[0])
        self.assertTrue(out.all() == b.all())

    def test_midpoints3_2(self):
        A = np.array([[1, 1, 1], [3, 3, 3], [5, 5, 5]])
        B = np.array([[2, 2, 2], [3, 3, 3], [5, 5, 5]])
        C = np.array([[3, 3, 3], [5, 5, 5], [7, 7, 7]])
        l1 = (A, B, C)
        a = tuple(midpoints(l1))
        out = np.array([[3.5, 3.5, 3.5], [2, 2, 2], [3, 3, 3], [5, 5, 5]])
        b = np.array(a[1])
        self.assertTrue(out.all() == b.all())

    def test_midpoints3_3(self):
        A = np.array([[1, 1, 1], [3, 3, 3], [5, 5, 5]])
        B = np.array([[2, 2, 2], [3, 3, 3], [5, 5, 5]])
        C = np.array([[3, 3, 3], [5, 5, 5], [7, 7, 7]])
        l1 = (A, B, C)
        a = tuple(midpoints(l1))
        out = np.array([[2, 2, 2], [3, 3, 3], [5, 5, 5],[4, 4, 4]])
        b = np.array(a[2])
        self.assertTrue(out.all() == b.all())

    def test_midpoints3_4(self):
        A = np.array([[1, 1, 1], [3, 3, 3], [5, 5, 5]])
        B = np.array([[2, 2, 2], [3, 3, 3], [5, 5, 5]])
        C = np.array([[3, 3, 3], [5, 5, 5], [7, 7, 7]])
        l1 = (A, B, C)
        a = tuple(midpoints(l1))
        out = np.array([[4, 4, 4],[3, 3, 3], [5, 5, 5], [7, 7, 7]])
        b = np.array(a[2])
        self.assertTrue(out.all() == b.all())