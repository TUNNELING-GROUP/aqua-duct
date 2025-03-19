 
import unittest
from unittest import TestCase
from aquaduct.utils.helpers import SmartRange
import numpy as np

class TestSmartRange(TestCase):
    def test_smartrange_array(self):
        itr = [2, 6, 6, 5, 4, 3, 7, 8, 9, 10, 11, 4, 4, 4, 4]
        S = SmartRange(itr)
        self.assertEqual(S.first_element(), 2)
        self.assertEqual(S.last_element(), 4)
        self.assertEqual(S.last_times(), 4)
        Sr = []
        for obj in S.raw:
            Sr.extend(obj.get())
        self.assertEqual(Sr, itr)
        self.assertEqual(list(S.get()), itr)
        S.append(5)
        itr.append(5)
        self.assertEqual(S.last_element(), 5)
        self.assertEqual(list(S.get()), itr)
        self.assertEqual(S.rev(), itr.reverse())
        
    def test_smartrange_ndarray(self):
        itr = np.array([2, 6, 6, 5, 4, 3, 7, 8, 9, 10, 11, 4, 4, 4, 4])
        S = SmartRange(itr)
        self.assertEqual(S.first_element(), 2)
        self.assertEqual(S.last_element(), 4)
        self.assertEqual(S.last_times(), 4)
        Sr = np.array([])
        for obj in S.raw:
            Sr = np.append(Sr, np.fromiter(obj.get(), int))
        for (elema, elemb) in zip(Sr, itr):
            self.assertEqual(elema, elemb)
        for (elema, elemb) in zip(list(S.get()), itr):
            self.assertEqual(elema, elemb)
        S.rev()
        for (elema, elemb) in zip(list(S.get()), np.flip(itr, 0)):
            self.assertEqual(elema, elemb)
        
    def test_smartrange_string(self):
        itr = 'aaaaaaaaaaaaaaabcdefgggggeeeeeeddijkfedcb'
        S = SmartRange(itr)
        self.assertEqual(S.first_element(), 'a')
        self.assertEqual(S.last_element(), 'b')
        self.assertEqual(S.last_times(), 1)
        itrs = ''
        for obj in S.raw:
            for elem in obj.get():
                itrs = itrs + elem
        self.assertEqual(itrs, itr)
        itrs = ''
        for obj in list(S.get()):
            itrs += obj
        self.assertEqual(itrs, itr)
             
        
           
if __name__ == '__main__':
    unittest.main()
        
