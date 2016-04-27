'''
Collection of helpers - functions and decorators.
'''

import numpy as np
from collections import Iterable
from functools import wraps
from os import close
from tempfile import mkstemp


########################################################################
#  aliens

def combine(seqin):
    '''returns a list of all combinations of argument sequences.
for example: combine((1,2),(3,4)) returns
[[1, 3], [1, 4], [2, 3], [2, 4]]'''

    # http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/302478/index_txt
    def rloop(seqin, listout, comb):
        '''recursive looping function'''
        if seqin:  # any more sequences to process?
            for item in seqin[0]:
                newcomb = comb + [item]  # add next item to current comb
                # call rloop w/ rem seqs, newcomb
                rloop(seqin[1:], listout, newcomb)
        else:  # processing last sequence
            listout.append(comb)  # comb finished, add to list

    listout = []  # listout initialization
    rloop(seqin, listout, [])  # start recursive process
    return listout


########################################################################

def lind(l, ind):
    """Indexes lists using lists of integers."""
    ll = []
    for i in ind:
        ll.append(l[i])
    return ll


class Auto:
    def __repr__(self):
        return "Auto"

    def __str__(self):
        return self.__repr__()


def create_tmpfile(ext=None):
    if ext is None:
        suffix = ''
    else:
        suffix = ".%s" % str(ext).lower()
    fd, name = mkstemp(suffix=suffix)
    close(fd)
    return name


def range2int(r, uniq=True):
    out = []
    for rr in r.split():
        if ':' in rr:
            if rr.count(':') == 1:
                r1, r2 = map(int, rr.split(':'))
                r3 = 1
            if rr.count(':') == 2:
                r1, r3, r2 = map(int, rr.split(':'))
            out.extend(range(r1, r2 + 1, r3))
        else:
            out.append(int(rr))
    if uniq:
        out = list(set(out))
        out.sort()
    return out


def int2range(l):
    '''
    Transforms list of integers in to a string of ranges.

    For example, a following list::

        [0,1,2,4,5,7,9]

    is transformed into::

        0:2 4:5 7 9

    :param list l: input list of int
    :return: string of ranges
    :rtype: str
    '''
    out = []
    l = list(set(l))
    l.sort()
    previous = None
    for e in l:
        if previous is None:
            previous = e
            out.append(e)
            out.append(':')
            out.append(None)
            continue
        else:
            if previous + 1 == e:
                out[-1] = e
                previous = e
                continue
            else:
                while out[-1] in [None, ':']:
                    out.pop(-1)
                out.append(' ')
                out.append(e)
                out.append(':')
                out.append(None)
                previous = e
                continue
    while out[-1] in [None, ':']:
        out.pop(-1)
    out = ''.join(map(str, out))
    return out


def is_iterable(l):
    '''
    Checks if provided obejct is iterable.
    Returns True is it is iterable, otherwise returns False.
    
    :param list l: input object
    
    :return: True if submited object is iterable otherwise returns False.
    :rtype: boolean

    .. warning::

        Current implementation cannot be used with generators!

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
        if isinstance(obj, Iterable):
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
        if isinstance(obj, Iterable):
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
        if isinstance(obj, Iterable):
            return np.matrix(list(obj)).A
        return np.matrix([obj]).A

    return patched


def list_blocks_to_slices(l):
    n = len(l)
    if n in [0, 1]:
        yield slice(None, None, None)
    if n > 1:
        prev = l[0]
        prev_nr = 0
        for nr, e in enumerate(l[1:]):
            if e == prev:
                continue
            yield slice(prev_nr, nr + 1, 1)
            prev = e
            prev_nr = nr + 1
        yield slice(prev_nr, nr + 2, 1)


@tupleify
def what2what(what, towhat):
    towhat = make_iterable(towhat)
    for nr, w in enumerate(make_iterable(what)):
        if w in towhat:
            yield nr


def make_iterable(something):
    if not is_iterable(something):
        return [something]
    return something
