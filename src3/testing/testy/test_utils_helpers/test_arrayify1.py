# -*- coding: utf-8 -*-
from unittest import TestCase
from aquaduct.utils.helpers import arrayify1
import numpy


class TestArrayify1(TestCase):
    def test_arrayify1(self):
        #chceck the type of output
        @arrayify1
        def do_sth_with_table(tables):
            tables = tables * 10
            return tables

        data = numpy.zeros((2, 3), float)
        the_stuff = do_sth_with_table(data)

        self.assertTrue(isinstance(the_stuff, numpy.ndarray))

    def test_size_of_matrix(self):
        # nupy arrays have size 1D=(x,)
        @arrayify1
        def do_sth_with_table(tables):
            for numbers in tables:
                yield '{0}'.format(numbers)

        tables = [[1, 1, 1, 1], [2, 2, 2, 2]]
        tables2 = [[1, 1], [2, 2]], [[3, 3], [4, 4]]
        res = (2,)
        res2 = (4,)
        obj = do_sth_with_table(tables)
        obj2 = do_sth_with_table(tables2)
        self.assertEqual(obj.shape, res)
        self.assertEqual(obj2.shape, res)


    def test_non_iterable(self):
        # non iterable output of internal function

        @arrayify1
        def do_sth_with_table(tables):
            output = tables * 10
            return (output)

        one_element = 10
        res = (1,)
        obj = do_sth_with_table(one_element)

        self.assertEqual(obj.shape, res)