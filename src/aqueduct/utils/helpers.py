'''
Collection of helpers - functions and decorators.
'''

from collections import Iterable
from functools import wraps

import numpy as np


def is_iterable(l):
    '''
    Checks if provided obejct is iterable.
    Returns True is it is iterable, otherwise returns False.
    
    :param list l: input object
    
    :return: True if submited object is iterable otherwise returns False.
    :rtype: boolean
    
    .. todo::
    
        Current implementation is primitive and HAVE TO be replaced.
    
    '''
    try:
        _ = (e for e in l)
        return True
    except TypeError:
        pass
    return False


def sortify(gen):
    '''
    Decorator to convert functions' outputs into a sorted list. If the output is iterable it is converted in to a list 
    of apropriate lenght. If the output is not iterable it is converted in to a list of lenght 1.
    
    Written on the basis of :func:`listify`.
    
    :returns: output of decorated function converted to a sorted list 
    :rtype: :py:class:`list`
    '''
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
    '''
    Decorator to convert functions' outputs into a list. If the output is iterable it is converted in to a list 
    of apropriate lenght. If the output is not iterable it is converted in to a list of lenght 1.
    
    This function was copied from:
    
    http://argandgahandapandpa.wordpress.com/2009/03/29/python-generator-to-list-decorator/
    
    and further improved by tljm@wp.pl.
    
    :returns: output of decorated function converted to a list 
    :rtype: :py:class:`list`
    '''
    @wraps(gen)
    def patched(*args, **kwargs):
        obj = gen(*args, **kwargs)
        if isinstance(obj,Iterable):
            return list(obj)
        return [obj]
    return patched

def tupleify(gen):
    '''
    Decorator to convert functions' outputs into a tuple. If the output is iterable it is converted in to a tuple 
    of apropriate lenght. If the output is not iterable it is converted in to a tuple of lenght 1.
    
    Written on the basis of :func:`listify`.
    
    :returns: output of decorated function converted to a tuple
    :rtype: :py:class:`tuple`
    '''
    @wraps(gen)
    def patched(*args, **kwargs):
        obj = gen(*args, **kwargs)
        if isinstance(obj,Iterable):
            return tuple(obj)
        return (obj,)
    return patched

def arrayify(gen):
    '''
    Decorator to convert functions' outputs into a 2D numpy array. If the output is iterable it is converted in to a a 2D numpy array
    of apropriate shape. If the output is not iterable it is converted in to a a 2D numpy array of shape 1x1.
    
    Written on the basis of :func:`listify`.
    
    :returns: output of decorated function converted to a 2D numpy array
    :rtype: :py:class:`numpy.ndarray`
    '''
    @wraps(gen)
    def patched(*args, **kwargs):
        obj = gen(*args, **kwargs)
        if isinstance(obj,Iterable):
            return np.matrix(list(obj)).A
        return np.matrix([obj]).A
    return patched
    
