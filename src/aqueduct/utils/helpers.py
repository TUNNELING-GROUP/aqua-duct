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
    '''
    This is an alien function. It is not extensively used.

    Directly taken form http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/302478/index_txt

    Returns a list of all combinations of argument sequences.
    For example, following call::

        combine(((1,2),(3,4)))

    gives following list of combinations::

        [[1, 3], [1, 4], [2, 3], [2, 4]]

    :param tuple sequin: Tuple of sequences to combine.
    :returns: All possible combinations of all input sequences.
    :rtype: list of lists
    '''
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
    """
    Indexes lists using lists of integers as identificators.
    For example::

        lind(['a','b','c','d','e'],[1,4,2])

    returns::

        ['b', 'e', 'c']

    :param list l: List to be indexed.
    :param list ind: Intiger indexes.
    :return: Reindexed list.
    :rtype: list

    """
    ll = []
    for i in ind:
        ll.append(l[i])
    return ll


class Auto:
    '''
    Auto type definition.
    The class is used as an alternative value for options (if particular option supports it).
    If options (or variables/parameters etc.) have value of :class:`Auto` it means that an automatic
    process for parametrization shoud be performed.

    For example, if the input parameter is set to :class:`Auto` it is supposed that its value is calculated
    on the basis of input data or other parapmeters.
    '''
    def __repr__(self):
        '''
        :return: String ``Auto``.
        :rtype: str
        '''
        return "Auto"

    def __str__(self):
        '''
        Calls :meth:`__repr__`.
        '''
        return self.__repr__()


def create_tmpfile(ext=None):
    '''
    Creates temporary file. File is created, closed and its file name is returned.

    .. note::

        It is responsability of the caller to delete the file.

    :param str ext: Optional extension of the file.
    :return: File name of created temporary file.
    :rtype: str
    '''
    if ext is None:
        suffix = ''
    else:
        suffix = ".%s" % str(ext).lower()
    fd, name = mkstemp(suffix=suffix)
    close(fd)
    return name


def range2int(r, uniq=True):
    '''
    Transforms a string range in to a list of integers.

    For example, a following string::

        '0:2 4:5 7 9'

    is transformed into::

        [0,1,2,4,5,7,9]

    :param str r: String of input range.
    :param bool uniq: Optional parameter, if set to `True` only unique and sorted integers are returned.
    :return: List of integers.
    :rtype: list of int
    '''
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
    Transforms a list of integers in to a string of ranges.

    For example, a following list::

        [0,1,2,4,5,7,9]

    is transformed into::

        0:2 4:5 7 9

    :param list l: input list of int
    :return: String of ranges.
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
    Checks if provided object is iterable.
    Returns True is it is iterable, otherwise returns False.

    :param list l: input object

    :return: True if submited object is iterable otherwise returns False.
    :rtype: bool

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

    :returns: Output of decorated function converted to a sorted list.
    :rtype: list
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

    :returns: Output of decorated function converted to a list.
    :rtype: list
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

    :returns: Output of decorated function converted to a tuple.
    :rtype: tuple
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

    :returns: Output of decorated function converted to a 2D numpy array.
    :rtype: numpy.ndarray
    '''

    @wraps(gen)
    def patched(*args, **kwargs):
        obj = gen(*args, **kwargs)
        if isinstance(obj, Iterable):
            return np.matrix(list(obj)).A
        return np.matrix([obj]).A

    return patched


def list_blocks_to_slices(l):
    '''
    Slices list in to block according to its elements identity. Resulting slices correspond to blocks of
    identical elements.

    :param list l: List of any objects.
    :return: Generator of slices.
    :rtype: generator
    '''
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
    '''
    what2what(what, towhat)
    :param list what: Input list for which indices of elements present in :attr:`towhat` are returned.
    :param list towhat: List of elements which input list is indexed to.
    :return: Indices of :attr:`what` list that are present in :attr:`towhat` list.
    :rtype: tuple
    '''
    towhat = make_iterable(towhat)
    for nr, w in enumerate(make_iterable(what)):
        if w in towhat:
            yield nr


def make_iterable(something):
    '''
    If input object is not interable returns it as one element list. Otherwise returns the object.

    :param object something: Input object.
    :return: Iterable object.
    :rtype: iterable or list
    '''
    if not is_iterable(something):
        return [something]
    return something
