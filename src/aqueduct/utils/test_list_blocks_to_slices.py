# coding: UTF-8
from unittest import TestCase
from helpers import list_blocks_to_slices
from helpers import listify

class TestList_blocks_to_slices(TestCase):
    #TODO jak przetestowac generator
    def test_how_to_slice(self):
        data=[1,1,1,2,2,2,'a','a','a']
        checks= listify(list_blocks_to_slices(data))
        #
        res=[[1,1,1],[2,2,2],['a','a','a']]
        self.assertEquals(checks,res)