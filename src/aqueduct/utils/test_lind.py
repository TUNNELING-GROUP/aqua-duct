from unittest import TestCase

def lind(l, ind):
    """
    Indexes lists using lists of integers as identificators.
    For example::

        lind(['a','b','c','d','e'],[1,4,2])

    returns::

        ['b', 'e', 'c']

    :param list l: List to be indexed.
    :param list ind: Intiger indexes.
    :return: Reindexed list.
    :rtype: list

    """
    ll = []
    for i in ind:
        ll.append(l[i])
    return ll

class TestLind(TestCase):
    pass


class TestLind(TestCase):
    def test_lind(self):
        testCase1= ["t","u","n","n","e","l]"]
        testCase2= [2,5,3]
        expected= ["u","l","n"]
        output=lind(testCase1,testCase2)
        if ( output != expected):
            Exception("test failed")

