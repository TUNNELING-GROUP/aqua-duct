from unittest import TestCase
from aqueduct.utils.helpers import tupleify

class TestTupleify(TestCase):
    def test_tupleify(self):
    # does it work?
        @tupleify
        def corvus_generator():
            corvus = ['jackdaw', 'magpie', 'raven', 'crow']
            for bird in corvus:
                yield '{0}'.format(bird)
        self.assertTrue(isinstance(tupleify(corvus_generator)(), tuple))

    def test_tupleify_string(self):
        # works even with string on input
        @tupleify
        def corvus_generator2():
            corvus = 'this string contains data about corvus'
            for bird in corvus:
                yield '{0}'.format(bird)

        self.assertTrue(isinstance(tupleify(corvus_generator2)(), tuple))

    def test_tupleify_tuple(self):
        # tupleify this tuple!
        @tupleify
        def corvus_generator3():
            corvus = ('jackdaw', 15, 'magpie', 3)
            for bird in corvus:
                yield '{0}'.format(bird)

        self.assertTrue(isinstance(tupleify(corvus_generator3)(), tuple))

    def test_tupleify_one_value(self):
        # numeric value can be tupleify
        @tupleify
        def corvus_generator4():
            corvus = 15, 17, 30
            for bird in corvus:
                yield '{0}'.format(bird)

        self.assertTrue(isinstance(tupleify(corvus_generator4)(), tuple))

    def test_correct_length_of_tuple(self):
        #if argument is list of tuples output is a tuple with number of elements equal to the number of
        #elements in list
        @tupleify
        def corvus_generator5(corvus):
            for bird in corvus:
                yield '{0}'.format(bird)
        black = ("raven", 'crow', 'jackdaw')
        blackandwhite = ("magpie")
        corvus = [black, blackandwhite]
        objectos=corvus_generator5(corvus)
        self.assertTrue(len(objectos)==2)
        self.assertTrue(isinstance(objectos,tuple))