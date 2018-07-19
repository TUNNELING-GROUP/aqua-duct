# -*- coding: utf-8 -*-
import types
from unittest import TestCase

from aquaduct.utils.helpers import list_blocks_to_slices


class TestList_blocks_to_slices(TestCase):

    def test_how_to_slice(self):
        data = [1, 1, 1, 2, 2, 2, 3, 3, 3]
        self.assertTrue(isinstance(list_blocks_to_slices(data), types.GeneratorType))

    def test_values(self):
        data = [1, 1, 1, 2, 2, 2, 3, 3, 3]
        checks = [data[sl] for sl in list_blocks_to_slices(data)]
        res = [[1, 1, 1], [2, 2, 2], [3, 3, 3]]
        self.assertEquals(checks, res)

    def test_values_other_types(self):
        # elements in input list can be different types, number of each elements in slices can differ
        dat = ["a", "a", 1, 1, 1, 2.3, 2.3]
        checks = [dat[sl] for sl in list_blocks_to_slices(dat)]
        res = [["a", "a"], [1, 1, 1], [2.3, 2.3]]
        self.assertEquals(checks, res)

    def test_slice_tuple(self):
        # if input is a tuple the function returns generator which creates tuples
        dat = (1, 1, 1, 'a', 'a', 'a', 'a')
        checks = [dat[sl] for sl in list_blocks_to_slices(dat)]
        res = [(1, 1, 1), ("a", "a", "a", "a")]
        self.assertEquals(checks, res)

    def test_slice_list_of_lists(self):
        # if input is a list of lists the output will be lists where internal of input lists will not be sliced
        dat = [[1, 2, 5], ['a', 'a', 'b', 'a']]
        checks = [dat[sl] for sl in list_blocks_to_slices(dat)]
        res = [[[1, 2, 5]], [['a', 'a', 'b', 'a']]]
        self.assertEquals(checks, res)
