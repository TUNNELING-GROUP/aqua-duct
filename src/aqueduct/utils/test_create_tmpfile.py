from unittest import TestCase
import shutil, tempfile
import os
from os import path
from helpers import create_tmpfile

# source: https://gist.github.com/odyniec/d4ea0959d4e0ba17a980

class TestCreate_tmpfile(TestCase):
    def setUp(self):
     # Create a temporary directory
        self.test_dir = tempfile.mkdtemp()

    def test_create_tmpfile(self):
        temporaryfile=create_tmpfile()
        tempproc=os.getpid()
        self.assertTrue(type(temporaryfile)==str)
        #how does the name of file is created?
        self.assertFalse(temporaryfile == tempproc)

    def test_chceck_directory(self):
        temporaryfile2 = create_tmpfile(ext='traj')
        self.assertEquals(temporaryfile2.endswith('.traj'), True)

