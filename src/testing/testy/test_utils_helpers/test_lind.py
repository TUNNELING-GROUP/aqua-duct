# -*- coding: utf-8 -*-
from unittest import TestCase
from aquaduct.utils.helpers import lind

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

