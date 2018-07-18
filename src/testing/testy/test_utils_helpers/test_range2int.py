# -*- coding: utf-8 -*-
from unittest import TestCase
from aquaduct.utils.helpers import range2int


class TestRange2int(TestCase):
    pass


class TestRange2int(TestCase):

    def test_letters(self):
        firstcase = range2int('0:2 4:5 7 9')
        firstres = [0, 1, 2, 4, 5, 7, 9]
        self.assertEquals(firstcase, firstres)

    def test_letters_below_zero(self):
        # if negative values are allowed
        fourthcase = range2int('-1:3 13:15')
        fourthres = [-1, 0, 1, 2, 3, 13, 14, 15]
        self.assertEquals(fourthcase, fourthres)

    def test_letters_cover_range(self):
        # if the ranges are covered in some part
        fivethcase = range2int('7:10 8:11')
        fivethres = [7, 8, 9, 10, 11]
        self.assertEquals(fivethcase, fivethres)

    def test_inverted_Range(self):
        # chceck if inverted range can break a call
        sixthcase = range2int("4:1 8:10")
        sixthres = [1, 2, 3, 4, 8, 9, 10]
        self.assertEquals(sixthres, sixthres)


class TestBadValue(TestCase):
    def test_letters_with_alpha(self):
        thirdcase = "a:10 15:10"
        self.assertRaises(ValueError, range2int, "a:10 15:10")

    def test_letters_and_commas(self):
        # no commas allowed
        self.assertRaises(ValueError, range2int, '1:2,4:5')


if __name__ == "__main__":
    unittest.main()
