# -*- coding: utf-8 -*-

import os

def get(filename):
        path = os.path.dirname(__file__)
        filename_with_path = os.path.join(path, filename)

        return os.path.abspath(filename_with_path)

