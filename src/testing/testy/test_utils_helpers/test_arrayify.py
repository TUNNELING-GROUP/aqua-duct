# -*- coding: utf-8 -*-

from unittest import TestCase
from aqueduct.utils.helpers import arrayify
import numpy

class TestArrayify(TestCase):

    def test_arrayify(self):
        @arrayify
        def do_sth_with_table(tables):
            tables=tables*10
            return tables

        data= numpy.zeros((2, 3), float)
        the_stuff=do_sth_with_table(data)

        self.assertTrue(isinstance(the_stuff, numpy.ndarray))

        #sprawdzić wymiary tablicy  shape

    def test_size_of_matrix(self):
        @arrayify
        def do_sth_with_table(tables):

            for numbers in tables:
                yield '{0}'.format(numbers)

        tables = [[1, 1, 1, 1], [2, 2, 2, 2]]
        tables2 = [[1, 1], [2, 2]], [[3, 3], [4, 4]]
        res=(1,2)
        res2=(2,2)
        obj=do_sth_with_table(tables)
        obj2=do_sth_with_table(tables2)
        self.assertEquals(obj.shape, res)
        self.assertEquals(obj2.shape, res)

    def test_size_other_input(self):
        #tuples and lists are good input
        @arrayify
        def do_sth_with_table(tables):
            for numbers in tables:
                yield '{0}'.format(numbers)
        tuple_in=((1,2,3),(4,5,6))
        res=(1,2)
        obj=do_sth_with_table(tuple_in)
        self.assertEquals(obj.shape, res)

    def test_non_iterable(self):
        # non iterable output of internal function
        @arrayify
        def do_sth_with_table(tables):
            output=tables*10
            return output

        one_element = 10
        res = (1, 1)
        obj = do_sth_with_table(one_element)

        self.assertEquals(obj.shape, res)