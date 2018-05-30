from aquaduct.utils.helpers import listify, SmartRange


################################################################################

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


################################################################################

@listify
def xor_full(a, b):
    ab = intersection_full(a, b)
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


################################################################################

def left_full(a, b):
    return intersection_full(a, xor_full(a, b))


def left_smartr(a, b):
    return intersection_smartr(a, xor_smartr(a, b))


def left_set(a, b):
    return intersection_set(a, xor_set(a, b))


################################################################################

def right_full(a, b):
    return intersection_full(b, xor_full(a, b))


def right_smartr(a, b):
    return intersection_smartr(b, xor_smartr(a, b))


def right_set(a, b):
    return intersection_set(b, xor_set(a, b))


################################################################################

intersection = intersection_set
left = left_set
right = right_set
