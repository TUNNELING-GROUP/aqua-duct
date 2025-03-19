import unittest
from aquaduct.utils.sets import presort

class TestPresort(unittest.TestCase):

    def test_presort(self):
        a = [1, 5]
        b = [3, 7]
        expected_s = [1, 5, 3, 7]
        expected_sa = [0, 2, 1, 3]
        s, sa = presort(a, b)
        self.assertEqual(s, expected_s)
        self.assertEqual(sa, expected_sa)

        a = [10, 20]
        b = [15, 25]
        expected_s = [10, 20, 15, 25]
        expected_sa = [0, 2, 1, 3]
        s, sa = presort(a, b)
        self.assertEqual(s, expected_s)
        self.assertEqual(sa, expected_sa)

        a = [5, 15]
        b = [10, 20]
        expected_s = [5, 15, 10, 20]
        expected_sa = [0, 2, 1, 3]
        s, sa = presort(a, b)
        self.assertEqual(s, expected_s)
        self.assertEqual(sa, expected_sa)

if __name__ == '__main__':
    unittest.main()
