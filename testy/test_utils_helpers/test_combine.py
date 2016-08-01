from unittest import TestCase
from aqueduct.utils.helpers import combine

class TestCombine(TestCase):
    pass


class TestCombineExist(TestCase):
    def test_combine(self):
        testing=((1,2),(3,4))
        case1=combine(testing)
        result1=[[1,3],[1,4],[2,3],[2,4]]
        self.assertTrue(case1==result1)

    def test_combine_failure(self):
        #does it take lists?
        testing = ([1, 2], [3, 4])
        case1_1 = combine(testing)
        result1_1 = [[1, 3], [1, 4], [2, 3], [2, 4]]
        self.assertTrue(case1_1== result1_1)

    def test_combine_example(self):
        testing2=((1, 2), (3, 4))
        case2 = combine(testing2)
        result2 = [[1, 2], [1, 4], [4, 3], [2, 4]]
        self.assertFalse(case2==result2)

    def test_combine_range(self):
        #the order of output is not so important
        testing2_2=((1,2),(3,4))
        case2_2=combine(testing2_2)
        result2_2=[[1,4],[1,3],[2,4],[2,3]]
        self.assertEquals(case2_2==result2_2)

    def test_combine_range(self):
        # function can take letters as argument
        testing2_2 = ((1, 'a'), (3, 4))
        case2_2 = combine(testing2_2)
        result2_2 = [[1, 3], [1, 4], ['a', 3], ['a', 4]]
        self.assertTrue(case2_2 == result2_2)