from aquaduct.utils.helpers import is_number

class SmartRangeFunction(object):
    __slots__ = "element times".split()

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

    def overlaps(self, srange):
        return (self.isin(srange.first_element()) or self.isin(srange.last_element()))

    def overlaps_mutual(self, srange):
        return self.overlaps(srange) or srange.overlaps(self)

    def contains(self, srange):
        # tests if srange of type SmartRange is in this range
        return self.isin(srange.first_element()) and self.isin(srange.last_element())


class SmartRangeEqual(SmartRangeFunction):
    __slots__ = "element times".split()
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
    __slots__ = "element times".split()
    type = 'i'

    def get(self):
        return (self.element + i for i in xrange(self.times))

    def rev(self):
        return SmartRangeDecrement(self.element + self.times - 1, self.times)

    def isin(self, element):
        return (element >= self.element) and (element <= self.element + self.times - 1)

    def last_element(self):
        return self.first_element() + self.times - 1


class SmartRangeDecrement(SmartRangeFunction):
    __slots__ = "element times".split()
    type = 'd'

    def get(self):
        return (self.element - i for i in xrange(self.times))

    def rev(self):
        return SmartRangeIncrement(self.element - self.times + 1, self.times)

    def isin(self, element):
        return (element <= self.element) and (element >= self.element - self.times + 1)

    def last_element(self):
        return self.first_element() - self.times + 1


class SmartRange(object):
    __slots__ = '__elements __len __min __max'.split()

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
        return '[%s]' % (','.join(map(str, self.__elements)))

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



if __name__ == "__main__":

    import sys
    from numbers import Number
    from collections import Set, Mapping, deque

    try:  # Python 2
        zero_depth_bases = (basestring, Number, xrange, bytearray)
        iteritems = 'iteritems'
    except NameError:  # Python 3
        zero_depth_bases = (str, bytes, Number, range, bytearray)
        iteritems = 'items'
    import re
    hidden = re.compile("^__")

    def getsize(obj_0):
        """Recursively iterate to sum size of object & members."""

        def inner(obj, _seen_ids=set()):
            obj_id = id(obj)
            if obj_id in _seen_ids:
                return 0
            _seen_ids.add(obj_id)
            size = sys.getsizeof(obj)
            if isinstance(obj, zero_depth_bases):
                pass  # bypass remaining control flow and return
            elif isinstance(obj, (tuple, list, Set, deque)):
                size += sum(inner(i) for i in obj)
            elif isinstance(obj, Mapping) or hasattr(obj, iteritems):
                size += sum(inner(k) + inner(v) for k, v in getattr(obj, iteritems)())
            # Check for custom object instances - may subclass above too
            if hasattr(obj, '__dict__'):
                size += inner(vars(obj))
            if hasattr(obj, '__slots__'):  # can have __slots__ with __dict__
                cn = "_%s__" % obj.__class__.__name__
                size += sum(inner(getattr(obj, s)) for s in obj.__slots__ if hasattr(obj, s))
                size += sum(inner(getattr(obj, re.sub(r"^__",cn,s))) for s in obj.__slots__ if hidden.match(s))
            return size

        return inner(obj_0)



    from aquaduct.utils import helpers

    l = 'fhjshfjkasdhfjkasdhfjkasdhfffffffffjkasdhfjkhsfruwerydbccvxcgdfffgffffgoooowoowwwwgfsdhgfhjasdgfhjasgdfsadgfgsdfhsdfsdfwerqwyruqyweuiyqweuiyq'*100

    s = helpers.SmartRange()
    s2 = SmartRange()

    print "old sr",getsize(s)
    print "new sr", getsize(s2)

    map(s.append,l)
    map(s2.append,l)

    del l

    print "old sr, l appended",getsize(s)
    print "new sr, l appended",getsize(s2)


    assert list(s.get()) == list(s2.get())

