# -*- coding: utf-8 -*-
from unittest import TestCase
from aqueduct.utils.helpers import is_iterable
from aqueduct.utils.helpers import Auto

class TestIs_iterable(TestCase):
    def test_is_iterable(self):
        pass


class TestIs_iterable(TestCase):
    def test_is_iterable(self):
        test1=["a","b","c","d"]
        self.assertTrue(is_iterable(test1))

    def test_is_iterable_num(self):
        test2 = [(1,2,3),(4,5,6)]
        self.assertTrue(is_iterable(test2))

    def test_is_iterable_not_list(self):
        test3 = ((1,2,3),(4,5,6))
        self.assertTrue(is_iterable(test3))

    def test_is_iterable_diff_type(self):
        test4 = [(1,2,3),"and now something completely different"]
        self.assertTrue(is_iterable(test4))

    def test_is_iterable_vector(self):
        test5 = (1,2,3,4,7,6,11)
        self.assertTrue(is_iterable(test5))

    def test_is_iterable_vector(self):
        test6 = (1, 2, 3, "a", 7, "b", 11)
        self.assertTrue(is_iterable(test6))

    def test_is_iterable_word(self):
        b="koty"
        test7 = b
        self.assertTrue(is_iterable(test7))

    def test_is_iterable_dictionary(self):
        test8 = {"a":"typ", "b":"zmiennej", "c": "jest", "d":"poprawny"}
        self.assertTrue(is_iterable(test8))

    def test_is_iterable_word(self):
        b =7
        test8 = b
        self.assertFalse(is_iterable(test8))

    def test_non_iterable(self):
        #class Auto is iterable
        obj=Auto()
        self.assertFalse(is_iterable(obj))