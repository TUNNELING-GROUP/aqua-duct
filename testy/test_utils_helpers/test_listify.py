# -*- coding: utf-8 -*-
from unittest import TestCase
from aqueduct.utils.helpers import listify



class TestListify(TestCase):
    pass

class TestingListyfy(TestCase):
    def test_if_list(self):
        #does it work?
        def corvus_generator():
            corvus = ['jackdaw', 'magpie', 'raven', 'crow']
            for bird in corvus:
                yield '{0}'.format(bird)

        self.assertTrue(isinstance(listify(corvus_generator)(), list))

    def test_isnot_list(self):
        # if function returns list it works well
        @listify
        def corvus_listificator(my_string):
            output=my_string.split(' ')
            return([output])
        corvus = 'jackdaw magpie raven crow'
        funk=corvus_listificator(corvus)
        self.assertTrue(isinstance(funk, list))

    def test_length(self):
        # if output of internal function is a list decorator returns just list of
        # output elements
        @listify
        def corvus_listificator(my_string):
            output = my_string.split(' ')
            return (output)
        corvus = 'jackdaw magpie raven crow'
        funk = corvus_listificator(corvus)
        self.assertTrue(len(funk)==4)