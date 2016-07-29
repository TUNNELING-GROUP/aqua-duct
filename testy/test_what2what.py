from unittest import TestCase
from aqueduct.utils.helpers import what2what

class TestWhat2what(TestCase):
    def test_what2what_output_type(self):
        list_what = [1, 2, 3]
        list_towhat = [1, 2, 2, 1, 3, 4]
        self.assertTrue(isinstance(what2what(list_what,list_towhat),tuple))

    def test_what2what(self):
        #returns empty tuple if element is not present
        #
        list_what=[50]
        list_what_2=[2,3]
        list_towhat=[1,2,2,1,3,4]
        check=what2what(list_what,list_towhat)
        check2=what2what(list_what_2, list_towhat)
        res=()
        res2=(0,1)
        self.assertEquals(check,res)
        self.assertEquals(check2,res2)

    def test_whatifoneisnotpresent(self):
        list_what_3 = [2, 50, 3]
        list_towhat=[1,2,2,3,4]
        check=what2what(list_what_3,list_towhat)
        res=(0,2)
        self.assertEquals(check,res)

class TestWat2what4diffTypes(TestCase):
    #function does not rise an exception if input are tuples
    def test_whatif_tuple(self):
        list_what=(1,2,3)
        list_towhat=(10,2,1,50,70)
        check=what2what(list_what,list_towhat)
        self.assertEqual(check,(0,1))

    def test_whatif_numeric(self):
        #if one list is just a number the output will be tuple (x,) where x is index from 'what'
        list_what = 4
        list_what_2=[2,50,1]
        list_towhat = (10, 2, 1, 4, 50, 70)
        list_towhat_2= 50
        check = what2what(list_what, list_towhat)
        check2=what2what(list_what_2,list_towhat_2)
        self.assertEqual(check, (0,))
        self.assertEquals(check2,(1,))

    def test_whatif_tuplevslist(self):
        #different types within input are allowed
        list_what=(1,3,4)
        list_what_2=[1,5]
        list_towhat=[40,10,5,1,100]
        list_towhat_2=(1,6,3,5)
        check=what2what(list_what,list_towhat)
        check2=what2what(list_what_2,list_towhat_2)
        res=(0,)
        res2=(0,1)
        self.assertEquals(check,res)
        self.assertEquals(check2,res2)

    def test_what_about_strins(self):
        #wow! it can deal with strings! :D
        list_what=("a","b","som")
        list_towhat="this string is iterable and it has 'b' inside, something "
        check=what2what(list_what,list_towhat)
        res=(0,1,2)
        self.assertEquals(check,res)
