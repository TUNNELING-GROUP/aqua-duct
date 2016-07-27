from unittest import TestCase
from helpers import tupleify

class TestTupleify(TestCase):
    def test_tupleify(self):
    # does it work?
        @tupleify
        def corvus_generator():
            corvus = ['jackdaw', 'magpie', 'raven', 'crow']
            for bird in corvus:
                yield '{0}'.format(bird)
        self.assertTrue(isinstance(tupleify(corvus_generator)(), tuple))
