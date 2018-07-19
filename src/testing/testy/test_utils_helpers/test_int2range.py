# -*- coding: utf-8 -*-
from unittest import TestCase

from aquaduct.utils.helpers import int2range


class TestInt2range(TestCase):
    def test_int2range(self):
        pass

    def test_range(self):
        firstcase = int2range([0, 1, 2, 4, 5, 7, 9])
        firstres = '0:2 4:5 7 9'
        self.assertEquals(firstcase, firstres)

    def test_below_zero(self):
        # if negative values are allowed
        fourthcase = int2range([-1, 0, 1, 2, 3, 13, 14, 15])
        fourthres = '-1:3 13:15'
        self.assertEquals(fourthcase, fourthres)


class TestBadValue(TestCase):

    def test_letters_with_alpha(self):
        arg1 = int2range([1, 2, 'a', 5])
        res = '1:2 5 a'
        self.assertEquals(arg1, res)
        # the function can order case with letter, letters are at the end of string

    def test_and_commas(self):
        # round squares allowed?
        self.assertTrue(int2range, (1, 2, 4, 5))
