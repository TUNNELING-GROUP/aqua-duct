import os


def get(filename):
    path = os.path.dirname(__file__)
    path = os.path.join(path, filename)

    return os.path.abspath(path)
