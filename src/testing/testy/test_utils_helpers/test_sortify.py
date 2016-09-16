# -*- coding: utf-8 -*-
from unittest import TestCase
from aqueduct.utils.helpers import sortify


def corvus_generator():
    corvus=['jackdaw','magpie','raven', 'crow']
    for bird in corvus:
        yield '{0}'.format(bird)

gen= corvus_generator()

class TestSortify(TestCase):
    def test_sortify(self):
        def corvus_generator():
            corvus = ['jackdaw', 'magpie', 'raven', 'crow']
            for bird in corvus:
                yield '{0}'.format(bird)

        self.assertTrue(isinstance(sortify(corvus_generator)(), list))

    def test_length(self):
        # if output of internal function is a string decorator returns just list of
        # output elements, if it is an empty string i returns one element
        @sortify
        def corvus_listificator(my_string):
            output = my_string.split(' ')
            return (output)
        corvus = 'jackdaw magpie raven crow'
        corvus2='jackdaw'
        corvus3= ''
        funk = corvus_listificator(corvus)
        funk2=corvus_listificator(corvus2)
        funk3=corvus_listificator(corvus3)

        self.assertTrue(len(funk)==4)
        self.assertTrue(len(funk2)==1)
        self.assertTrue(len(funk3)==1)

    def test_sort(self):
        # check if elements are sorted
        @sortify
        def corvus_list(my_string):
            output = my_string.split(' ')
            return (output)

        corvus = 'jackdaw magpie raven crow'
        funk = corvus_list(corvus)
        expect = ['crow', 'jackdaw', 'magpie', 'raven']
        self.assertTrue(funk == expect)

    def test_numbers(self):
        @sortify
        def number_list_from_string(my_String):
            output=my_String.split(' ')
            return (output)
        listing_numbers='1 5 7 3 4 -1 0 0 5'
        funky=number_list_from_string(listing_numbers)
        expect=['-1','0','0','1','3','4','5','5','7']
        self.assertTrue(funky==expect)