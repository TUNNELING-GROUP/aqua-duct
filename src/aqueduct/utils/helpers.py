'''
Collection of helpers - functions and decorators.

Decorators:

    .. autosummary::
    
        listify
        sortify
        tupleify
        arraify
        
Functions:

    .. autosummary::
    
        is_iterable

'''

from collections import Iterable
from functools import wraps

import numpy as np





def is_iterable(l):
    '''
    Checks if provided obejct is iterable.
    Returns True is it is iterable, otherwise returns False.
    
        :param str sender
    
    '''
    try:
        _ = (e for e in l)
        return True
    except TypeError:
        pass
    return False


def sortify(gen):
    # http://argandgahandapandpa.wordpress.com/2009/03/29/python-generator-to-list-decorator/
    # improved by tljm: for non iterable objects it returns list of one element: the object
    @wraps(gen)
    def patched(*args, **kwargs):
        obj = gen(*args, **kwargs)
        if is_iterable(obj):
            obj = list(obj)
            obj.sort()
            return obj
        return [obj]

    return patched

def listify(gen):
    '''Convert a generator into a function which returns a list
    #http://argandgahandapandpa.wordpress.com/2009/03/29/python-generator-to-list-decorator/
    # improved by TM: for non iterable objects it returns list of one element: the object
    '''
    @wraps(gen)
    def patched(*args, **kwargs):
        obj = gen(*args, **kwargs)
        if isinstance(obj,Iterable):
            return list(obj)
        return [obj]
    return patched

def tupleify(gen):
    @wraps(gen)
    def patched(*args, **kwargs):
        obj = gen(*args, **kwargs)
        if isinstance(obj,Iterable):
            return tuple(obj)
        return (obj,)
    return patched

def arrayify(gen):
    @wraps(gen)
    def patched(*args, **kwargs):
        obj = gen(*args, **kwargs)
        if isinstance(obj,Iterable):
            return np.matrix(list(obj)).A
        return np.matrix([obj]).A
    return patched
    