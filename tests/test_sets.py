import unittest

from aquaduct.utils.sets import (
    presort, intersection_simple, intersection_full, intersection_smartr, intersection_set,
    glue, glue_simple, xor_full, xor_smartr, xor_set, xor_simple,
    left_full, left_smartr, left_set, left_simple,
    right_full, right_smartr, right_set, right_simple
)

class TestSets(unittest.TestCase):

    def test_presort(self):
        self.assertEqual(presort([1, 3], [2, 4]), ([1, 3, 2, 4], [0, 2, 1, 3]))

    def test_intersection_simple(self):
        self.assertEqual(intersection_simple([1, 3], [2, 4]), [2, 3])
        self.assertEqual(intersection_simple([1, 2], [3, 4]), [])

    def test_intersection_full(self):
        self.assertEqual(intersection_full([1, 2, 3], [2, 3, 4]), [2, 3])

    def test_intersection_smartr(self):
        self.assertEqual(intersection_smartr([1, 2, 3], [2, 3, 4]), [2, 3])

    def test_intersection_set(self):
        self.assertEqual(intersection_set([1, 2, 3], [2, 3, 4]), [2, 3])

    def test_glue(self):
        self.assertEqual(glue([1, 2, 3], [3, 4, 5]), [1, 2, 3, 4, 5])
        self.assertEqual(glue([1, 2, 3], [4, 5, 6]), [])

    def test_glue_simple(self):
        self.assertEqual(glue_simple([1, 2, 3], [3, 4, 5]), [1, 2, 3, 4, 5])
        self.assertEqual(glue_simple([1, 2, 3], [4, 5, 6]), [])

    def test_xor_full(self):
        self.assertEqual(xor_full([1, 2, 3], [2, 3, 4]), [1, 4])

    def test_xor_smartr(self):
        self.assertEqual(xor_smartr([1, 2, 3], [2, 3, 4]), [1, 4])

    def test_xor_set(self):
        self.assertEqual(xor_set([1, 2, 3], [2, 3, 4]), [1, 4])

    def test_xor_simple(self):
        self.assertEqual(xor_simple([1, 2, 3], [2, 3, 4]), [1, 4])

    def test_left_full(self):
        self.assertEqual(left_full([1, 2, 3], [2, 3, 4]), [1])

    def test_left_smartr(self):
        self.assertEqual(left_smartr([1, 2, 3], [2, 3, 4]), [1])

    def test_left_set(self):
        self.assertEqual(left_set([1, 2, 3], [2, 3, 4]), [1])

    def test_left_simple(self):
        self.assertEqual(left_simple([1, 2, 3], [2, 3, 4]), [1])

    def test_right_full(self):
        self.assertEqual(right_full([1, 2, 3], [2, 3, 4]), [4])

    def test_right_smartr(self):
        self.assertEqual(right_smartr([1, 2, 3], [2, 3, 4]), [4])

    def test_right_set(self):
        self.assertEqual(right_set([1, 2, 3], [2, 3, 4]), [4])

    def test_right_simple(self):
        self.assertEqual(right_simple([1, 2, 3], [2, 3, 4]), [4])

if __name__ == '__main__':
    unittest.main()
