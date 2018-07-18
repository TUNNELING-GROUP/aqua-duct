# -*- coding: utf-8 -*-
from unittest import TestCase
from aquaduct.utils.helpers import make_iterable


class TestMake_iterable(TestCase):
    def test_type(self):
        iterable = (('tuple'), ('of'), ('tuples'))
        iterable2 = [('list'), ('of'), ['tuples and list']]
        iterable3 = "it is said that strings are iterable"
        noniterable = 1
        check1 = make_iterable(iterable)
        check2 = make_iterable(iterable2)
        check3 = make_iterable(iterable3)
        chceck4 = make_iterable(noniterable)

        self.assertEquals(check1, iterable)
        self.assertEquals(check2, iterable2)
        self.assertEquals(check3, iterable3)
        self.assertEquals(chceck4, [noniterable])

    def test_empty(self):
        # function will throw an exception if something is empty
        check = ()
        self.assertRaises(TypeError, make_iterable(check))
