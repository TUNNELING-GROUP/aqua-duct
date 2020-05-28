from aquaduct.utils.helpers import listify, SmartRange
import numpy as np


################################################################################

def presort(a, b):
    s = [a[0], a[-1], b[0], b[-1]]
    return s, np.argsort(s).tolist()


def intersection_simple(a, b):
    s, sa = presort(a, b)
    ss = sorted(sa[1:3])
    if sa in [[0, 1, 2, 3], [2, 3, 0, 1], [0, 2, 1, 3], [2, 0, 3, 1], [0, 2, 3, 1], [2, 0, 1, 3]]:
        if sa[1:3] in [[1, 2], [3, 0]]:
            if s[1] == s[2]:
                return list(range(s[sa[1]], s[sa[2]] + 1))
        else:
            return list(range(s[sa[1]], s[sa[2]] + 1))
    return []


def intersection_full(a, b):
    return [aa for aa in a if aa in b]


def intersection_smartr(a, b):
    b_ = SmartRange(b)
    return [aa for aa in a if b_.isin(aa)]


def intersection_set(a, b):
    return sorted(set(a).intersection(b))


################################################################################

def glue(a, b):
    if a[-1] >= b[0]:
        return sorted(set(a + b))
    return []


def glue_simple(a, b):
    if a[-1] >= b[0]:
        if a[0] > b[-1]:
            return list(range(min(a[0], b[0]), b[-1] + 1)) + list(range(a[0], max(a[-1], b[-1]) + 1))
        return list(range(min(a[0], b[0]), max(a[-1], b[-1]) + 1))
    return []


################################################################################

@listify
def xor_full(a, b):
    ab = set(intersection_full(a, b))
    for e in glue(a, b):
        if e not in ab:
            yield e


@listify
def xor_smartr(a, b):
    ab = SmartRange(intersection_smartr(a, b))
    for e in glue(a, b):
        if not ab.isin(e):
            yield e


def xor_set(a, b):
    ab = set(intersection_set(a, b))
    return [e for e in glue(a, b) if e not in ab]


def xor_simple(a, b):
    lohi = intersection_simple(a, b)
    if len(lohi):
        lo, hi = lohi[0], lohi[-1]
        return [e for e in glue_simple(a, b) if e < lo or e > hi]
    return glue_simple(a, b)


################################################################################

def left_full(a, b):
    return intersection_full(a, xor_full(a, b))


def left_smartr(a, b):
    return intersection_smartr(a, xor_smartr(a, b))


def left_set(a, b):
    return intersection_set(a, xor_set(a, b))


def left_simple(a, b):
    # return intersection_simple(a, xor_simple(a, b))
    # intersection simple is too simple!
    return intersection_set(a, xor_simple(a, b))


################################################################################

def right_full(a, b):
    return intersection_full(b, xor_full(a, b))


def right_smartr(a, b):
    return intersection_smartr(b, xor_smartr(a, b))


def right_set(a, b):
    return intersection_set(b, xor_set(a, b))


def right_simple(a, b):
    # intersection simple is too simple!
    return intersection_set(b, xor_simple(a, b))


################################################################################

intersection = intersection_set
left = left_simple
right = right_simple
