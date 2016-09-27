from unittest import TestCase
from aqueduct.geom.traces import LinearizeOneWay

class TestLinearizeOneWay(TestCase):
    def test_Linearize_One_Way(self):
        test_case=((2,3,4),(2,3,6),(2,3,8))
        self.assertTrue(True==test_case.here())

