# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2016-2018  Tomasz Magdziarz, Alicja Płuciennik, Michał Stolarczyk <info@aquaduct.pl>
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
from functools import partial,total_ordering
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

def robust_and(a,b):
    if a is None:
        return bool(b)
    if b is None:
        return bool(a)
    return a and b

def robust_or(a,b):
    if a is None:
        return bool(b)
    if b is None:
        return bool(a)
    return a or b


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

def noaction(gen):

    @wraps(gen)
    def patched(*args, **kwargs):
        return gen(*args, **kwargs)

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

class arrayify(object):

    def __init__(self,shape=None):
        self.shape = shape

    def __call__(self,gen):
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
                result = np.matrix(list(obj),dtype=defaults.float_default).A
            else:
                result = np.matrix([obj],dtype=defaults.float_default).A
            if self.shape is None: return result
            new_shape = []
            for ds,s in zip(self.shape,result.shape):
                if ds is None:
                    if result.size:
                        new_shape.append(s)
                    else:
                        new_shape.append(0)
                else:
                    new_shape.append(ds)
            result.shape = tuple(new_shape)
            return result

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

def iterate_or_die(something,times=None):
    if is_iterable(something):
        return something
    return (something for dummy in xrange(times))


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

#@total_ordering
class SmartRangeFunction(object):
    def __init__(self, element, times):
        self.element = element
        self.times = times

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return "%s(%r,%d)" % (self.__class__.__name__, self.element, self.times)

    def __len__(self):
        return self.times
    
    def get(self):
        raise NotImplementedError('This method should be implemented in a child class.')

    def rev(self):
        raise NotImplementedError('This method should be implemented in a child class.')

    def isin(self, element):
        raise NotImplementedError('This method should be implemented in a child class.')

    def first_element(self):
        return self.element

    def last_element(self):
        # this is suboptimal, implement it in child class
        return self.get()[-1]

    def overlaps(self,srange):
        return (self.isin(srange.first_element()) or self.isin(srange.last_element()))

    def overlaps_mutual(self,srange):
        return self.overlaps(srange) or srange.overlaps(self)

    def contains(self,srange):
        # tests if srange of type SmartRange is in this range
        return self.isin(srange.first_element()) and self.isin(srange.last_element()) 



class SmartRangeEqual(SmartRangeFunction):

    type = 'e'

    def get(self):
        return [self.element] * self.times

    def rev(self):
        return self

    def isin(self, element):
        return element == self.element

    def last_element(self):
        return self.first_element()


class SmartRangeIncrement(SmartRangeFunction):

    type = 'i'

    def get(self):
        return (self.element + i for i in xrange(self.times))

    def rev(self):
        return SmartRangeDecrement(self.element + self.times - 1, self.times)

    def isin(self, element):
        return (element >= self.element) and (element <= self.element + self.times - 1)

    def last_element(self):
        return self.first_element()+self.times-1


class SmartRangeDecrement(SmartRangeFunction):

    type = 'd'

    def get(self):
        return (self.element - i for i in xrange(self.times))

    def rev(self):
        return SmartRangeIncrement(self.element - self.times + 1, self.times)

    def isin(self, element):
        return (element <= self.element) and (element >= self.element - self.times + 1)

    def last_element(self):
        return self.first_element()-self.times+1

class SmartRange(object):
    def __init__(self, iterable=None):
        self.__elements = []
        self.__len = 0
        self.__min = None
        self.__max = None

        if iterable is not None:
            map(self.append, iterable)

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return '[%s]' % (','.join(map(str,self.__elements)))


    def first_element(self):
        if len(self.__elements) == 0:
            return None
        element = self.__elements[0]
        if isinstance(element, SmartRangeFunction):
            return element.element
        return element

    def last_element(self):
        if len(self.__elements) == 0:
            return None
        element = self.__elements[-1]
        if isinstance(element, SmartRangeFunction):
            return element.element
        return element

    def last_times(self):
        if len(self.__elements) == 0:
            return 0
        element = self.__elements[-1]
        if isinstance(element, SmartRangeFunction):
            return element.times
        return 1

    @property
    @listify
    def raw(self):
        for element in self.__elements:
            if not isinstance(element, SmartRangeFunction):
                yield SmartRangeEqual(element, 1)
            else:
                yield element

    def append(self, element):
        assert not isinstance(element, SmartRangeFunction)
        if len(self.__elements) == 0:
            self.__elements.append(element)
            self.__min = element
            self.__max = element
        else:
            if element == self.last_element():
                if isinstance(self.__elements[-1], SmartRangeEqual) or (
                        not isinstance(self.__elements[-1], SmartRangeFunction)):
                    self.__elements[-1] = SmartRangeEqual(element, self.last_times() + 1)
                else:
                    self.__elements.append(element)
            else:
                if not is_number(element):
                    self.__elements.append(element)
                else:
                    if element - self.last_times() == self.last_element():
                        if isinstance(self.__elements[-1], SmartRangeIncrement) or (
                                not isinstance(self.__elements[-1], SmartRangeFunction)):
                            self.__elements[-1] = SmartRangeIncrement(self.last_element(), self.last_times() + 1)
                        else:
                            self.__elements.append(element)
                    elif element + self.last_times() == self.last_element():
                        if isinstance(self.__elements[-1], SmartRangeDecrement) or (
                                not isinstance(self.__elements[-1], SmartRangeFunction)):
                            self.__elements[-1] = SmartRangeDecrement(self.last_element(), self.last_times() + 1)
                        else:
                            self.__elements.append(element)
                    else:
                        self.__elements.append(element)
            if element > self.__max:
                self.__max = element
            if element < self.__min:
                self.__min = element
        self.__len += 1

    def get(self):
        for element in self.__elements:
            if not isinstance(element, SmartRangeFunction):
                yield element
            else:
                for e in element.get():
                    yield e

    def rev(self):
        elements = []
        for e in self.__elements[::-1]:
            if isinstance(e, SmartRangeFunction):
                elements.append(e.rev())
            else:
                elements.append(e)
        self.__elements = elements

    def __len__(self):
        return self.__len

    def __iter__(self):
        return self.get()

    def min(self):
        return self.__min

    def max(self):
        return self.__max

    def isin(self, element):
        for block in self.raw:
            if block.isin(element):
                return True
        return False


