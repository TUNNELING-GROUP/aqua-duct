# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2016-2017  Tomasz Magdziarz, Alicja Płuciennik, Michał Stolarczyk <info@aquaduct.pl>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Collection of helpers - functions and decorators.
"""

import numpy as np
from collections import Iterable
from functools import wraps
from os import close
from tempfile import mkstemp

from aquaduct.utils.maths import defaults


########################################################################
#  aliens

def combine(seqin):
    """
    This is an alien function. It is not extensively used.

    Directly taken form http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/302478/index_txt

    Returns a list of all combinations of argument sequences.
    For example, following call::

        combine(((1,2),(3,4)))

    gives following list of combinations::

        [[1, 3], [1, 4], [2, 3], [2, 4]]

    :param tuple seqin: Tuple of sequences to combine.
    :returns: All possible combinations of all input sequences.
    :rtype: list of lists
    """

    # http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/302478/index_txt
    def rloop(seqin, listout, comb):
        """recursive looping function"""
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


def are_rows_uniq(some_array):
    #http://stackoverflow.com/questions/16970982/find-unique-rows-in-numpy-array
    ca = np.ascontiguousarray(some_array).view(np.dtype((np.void, some_array.dtype.itemsize * some_array.shape[1])))
    return np.unique(ca).shape[0] == ca.shape[0]


########################################################################

def is_number(s):
    # http://pythoncentral.org/how-to-check-if-a-string-is-a-number-in-python-including-unicode/
    if isinstance(s, bool):
        return False
    try:
        float(s)
        return True
    except:
        pass
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except:
        pass
    return False


def lind(l, ind):
    """
    Indexes lists using lists of integers as identificators.
    For example::

        lind(['a','b','c','d','e'],[1,4,2])

    returns::

        ['b', 'e', 'c']

    :param list l: List to be indexed.
    :param list ind: Integer indexes.
    :return: Reindexed list.
    :rtype: list

    """
    ll = []
    for i in ind:
        ll.append(l[i])
    return ll


class Auto:
    """
    Auto type definition.
    The class is used as an alternative value for options (if particular option supports it).
    If options (or variables/parameters etc.) have value of :class:`Auto` it means that an automatic
    process for parametrization should be performed.

    For example, if the input parameter is set to :class:`Auto` it is supposed that its value is calculated
    on the basis of input data or other parameters.
    """

    def __repr__(self):
        """
        :return: String ``Auto``.
        :rtype: str
        """
        return "Auto"

    def __str__(self):
        """
        Calls :meth:`__repr__`.
        """
        return self.__repr__()


def create_tmpfile(ext=None):
    """
    Creates temporary file. File is created, closed and its file name is returned.

    .. note::

        It is responsibility of the caller to delete the file.

    :param str ext: Optional extension of the file.
    :return: File name of created temporary file.
    :rtype: str
    """
    if ext is None:
        suffix = ''
    else:
        suffix = ".%s" % str(ext).lower()
    fd, name = mkstemp(suffix=suffix)
    close(fd)
    return name


def range2int(r, uniq=True):
    """
    Transforms a string range in to a list of integers (with added missing elements from given ranges).

    For example, a following string::

        '0:2 4:5 7 9'

    is transformed into::

        [0,1,2,4,5,7,9]

    :param str r: String of input range.
    :param bool uniq: Optional parameter, if set to `True` only unique and sorted integers are returned.
    :return: List of integers.
    :rtype: list of int
    """
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
    # the function can order case with letter, letters are at the end of string
    """
    Transforms a list of integers in to a string of ranges.

    For example, a following list::

        [0,1,2,4,5,7,9]

    is transformed into::

        0:2 4:5 7 9

    :param list l: input list of int
    :return: String of ranges.
    :rtype: str
    """
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
    """
    Checks if provided object is iterable.
    Returns True is it is iterable, otherwise returns False.

    :param list l: input object

    :return: True if submitted object is iterable otherwise returns False.
    :rtype: bool

    .. warning::

        Current implementation cannot be used with generators!

    .. todo::

        Current implementation is primitive and HAVE TO be replaced.

    """
    try:
        _ = (e for e in l)
        return True
    except TypeError:
        pass
    return False


def sortify(gen):
    """
    Decorator to convert functions' outputs into a sorted list. If the output is iterable it is converted in to a list
    of appropriate length. If the output is not iterable it is converted in to a list of length 1.

    Written on the basis of :func:`listify`.

    :returns: Output of decorated function converted to a sorted list.
    :rtype: list
    """

    @wraps(gen)
    def patched(*args, **kwargs):
        obj = gen(*args, **kwargs)
        if is_iterable(obj):
            obj = list(obj)
            obj.sort()
            return obj
        return [obj]

    return patched


def uniqify(gen):
    """
    Decorator to convert functions' outputs into a sorted list of unique objects. If the output is iterable it is
    converted in to a list of appropriate length. If the output is not iterable it is converted in to a list of length 1.

    Written on the basis of :func:`listify`.

    :returns: Output of decorated function converted to a sorted list of unique objects.
    :rtype: list
    """

    @wraps(gen)
    def patched(*args, **kwargs):
        obj = gen(*args, **kwargs)
        if is_iterable(obj):
            obj = list(set(obj))
            obj.sort()
            return obj
        return [obj]

    return patched


def listify(gen):
    """
    Decorator to convert functions' outputs into a list. If the output is iterable it is converted in to a list
    of appropriate length. If the output is not iterable it is converted in to a list of length 1.

    This function was copied from:

    http://argandgahandapandpa.wordpress.com/2009/03/29/python-generator-to-list-decorator/

    and further improved by tljm@wp.pl.

    :returns: Output of decorated function converted to a list.
    :rtype: list
    """

    @wraps(gen)
    def patched(*args, **kwargs):
        obj = gen(*args, **kwargs)
        if isinstance(obj, Iterable):
            return list(obj)
        return [obj]

    return patched


def tupleify(gen):
    """
    Decorator to convert functions' outputs into a tuple. If the output is iterable it is converted in to a tuple
    of apropriate length. If the output is not iterable it is converted in to a tuple of length 1.

    Written on the basis of :func:`listify`.

    :returns: Output of decorated function converted to a tuple.
    :rtype: tuple
    """

    @wraps(gen)
    def patched(*args, **kwargs):
        obj = gen(*args, **kwargs)
        if isinstance(obj, Iterable):
            return tuple(obj)
        return (obj,)

    return patched


def arrayify(gen):
    """
    Decorator to convert functions' outputs into a 2D numpy array. If the output is iterable it is converted in to a 2D numpy array
    of appropriate shape. If the output is not iterable it is converted in to a 2D numpy array of shape 1x1.

    Written on the basis of :func:`listify`.

    :returns: Output of decorated function converted to a 2D numpy array.
    :rtype: numpy.ndarray
    """

    @wraps(gen)
    def patched(*args, **kwargs):
        obj = gen(*args, **kwargs)
        if isinstance(obj, Iterable):
            return np.matrix(list(obj),dtype=defaults.float_default).A
        return np.matrix([obj],dtype=defaults.float_default).A

    return patched


def arrayify1(gen):
    """
    Decorator to convert functions' outputs into a 1D numpy array. If the output is iterable it is converted in to a 2D numpy array
    of appropriate shape. If the output is not iterable it is converted in to a 2D numpy array of shape 1x1.

    Written on the basis of :func:`listify`.

    :returns: Output of decorated function converted to a 1D numpy array.
    :rtype: numpy.ndarray
    """

    @wraps(gen)
    def patched(*args, **kwargs):
        obj = gen(*args, **kwargs)
        if isinstance(obj, Iterable):
            return np.matrix(list(obj),dtype=defaults.float_default).A1
        return np.matrix([obj],dtype=defaults.float_default).A1

    return patched


def list_blocks_to_slices(l):
    """
    Slices list in to block according to its elements identity. Resulting slices correspond to blocks of
    identical elements.

    :param list l: List of any objects.
    :return: Generator of slices.
    :rtype: generator
    """
    # TODO: poprawic opis
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

def split_list(l,s):
    # l is list
    # s is element to split
    if s not in l:
        yield l
    else:
        n = l.count(s)
        i = -1
        b = 0
        while n:
            i = l.index(s,i+1)
            yield l[b:i]
            b = i+1
            n -= 1
        yield l[b:]



@tupleify
def what2what(what, towhat):
    """
    what2what(what, towhat)
    This function search if elements of the one list (:attr: 'what') are present in the other list (:attr: 'towhat') and returns indices of elements form :attr:'what' list as a tuple.
    If elements from the first list are not present in the second list the tuple is empty.
    :param list what: Input list for which indices of elements present in :attr:`towhat` are returned.
    :param list towhat: List of elements which input list is indexed to.
    :return: Indices of :attr:`what` list that are present in :attr:`towhat` list.
    :rtype: tuple
    """
    # todo poprawic opis
    towhat = make_iterable(towhat)
    for nr, w in enumerate(make_iterable(what)):
        if w in towhat:
            yield nr


def make_iterable(something):
    """
    If input object is not iterable returns it as one element list. Otherwise returns the object.

    :param object something: Input object.
    :return: Iterable object.
    :rtype: iterable or list
    """
    if not is_iterable(something):
        return [something]
    return something


def strech_zip(*args):
    ns = map(float, map(len, args))
    N = int(max(ns))
    for n in range(N):
        yield tuple([args[nr][int(cN / N * n)] for nr, cN in enumerate(ns)])


def compress_zip(*args):
    ns = map(float, map(len, args))
    N = int(min(ns))
    position = [0.] * len(args)
    for n in range(N):
        this_yield = []
        next_position = [float(len(a)) / N + p for a, p in zip(args, position)]
        for a, p, np in zip(args, position, next_position):
            if n + 1 == N:
                this_yield.append(a[int(p):])
            else:
                this_yield.append(a[int(p):int(np)])
        yield tuple(this_yield)
        position = next_position


def zip_zip(*args, **kwargs):
    if 'N' in kwargs.keys():
        N = kwargs['N']
    else:
        N = int(min(map(float, map(len, args))))
    position = [0.] * len(args)
    for n in range(N):
        this_yield = []
        next_position = [float(len(a)) / N + p for a, p in zip(args, position)]
        for a, p, np in zip(args, position, next_position):
            ip = int(p)
            inp = int(np)
            if n + 1 == N:
                this_yield.append(a[ip:])
            else:
                if ip == inp:
                    inp += 1
                this_yield.append(a[ip:inp])
        yield tuple(this_yield)
        position = next_position


def xzip_xzip(*args, **kwargs):
    if 'N' in kwargs.keys():
        N = kwargs['N']
    else:
        N = int(min(map(float, args)))
    position = [0.] * len(args)

    for n in xrange(N):
        this_yield = []
        # next_position = [float(a) / N + p for a, p in zip(args, position)]
        next_position_ = (float(args[i]) / N + position[i] for i in xrange(len(args)))
        next_position = []
        for i in xrange(len(args)):
            a = args[i]
            ip = int(position[i])
            next_position.append(next_position_.next())
            inp = int(next_position[-1])
            if n + 1 == N:
                this_yield.append(slice(ip, None))
            else:
                if ip == inp:
                    inp += 1
                this_yield.append(slice(ip, inp))
        yield tuple(this_yield)
        position = next_position

def concatenate(*args):
    '''
    Concatenates input iterable arguments in to one generator.
    '''
    for a in args:
        for e in a:
            yield e

class Bunch(object):
    """
    http://code.activestate.com/recipes/52308
    foo=Bunch(a=1,b=2)
    """
    def __init__(self, **kwds):
        self.__dict__.update(kwds)
