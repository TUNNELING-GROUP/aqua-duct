# -*- coding: utf-8 -*-



import numpy as np

import MDAnalysis as mda

from aqueduct.geom.convexhull import ConvexHull
from aqueduct.utils.helpers import int2range


class Selection(object):
    '''
    def __init__(self,selection,selection_string=None):

        self.selection_object = selection
        self.selection_string = selection_string
    '''

    def center_of_mass(self):
        # should return numpy (3,) array
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def iterate_over_residues(self):
        # should iterate over residues in the selection returning object of the same type
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def unique_resids(self):
        # should return array of resids
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def unique_resids_number(self):
        return len(self.unique_resids(ikwid=True))

    def atom_positions(self):
        # should return numpy (x,3) array
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def center_of_mass_of_residues(self):
        # should resturn list of lists or generator of center of masses
        return (R.center_of_mass().tolist() for R in self.iterate_over_residues())

    def get_convexhull_of_atom_positions(self):
        # should return modified ConvexHull object
        return ConvexHull(self.atom_positions())

    def uniquify(self):
        # should change selection to unique atoms only
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def __add__(self, other):
        raise NotImplementedError("This is abstract class. Missing implementation in a child class.")

    def first_resid(self):
        return self.resids.tolist()[0]


# TODO: decide if methods should be properties or not

class SelectionMDA(mda.core.AtomGroup.AtomGroup,
                   Selection):
    def iterate_over_residues(self):
        return (self.__class__(R) for R in self.residues)

    def unique_resids(self, ikwid=False):
        # TODO: do something with this method!
        assert ikwid, "This causes bugs! Avoid this method or take special care in using its results. If you want to use it pass additional variable ikwid = True."
        return np.unique(self.resids)

    def atom_positions(self):
        return self.positions

    def __add__(self, other):
        return SelectionMDA(self._atoms + other._atoms)

    def uniquify(self):
        self.__init__(mda.core.AtomGroup.AtomGroup(set(self._atoms)))


class CompactSelectionMDA(object):
    def __init__(self, sMDA):

        self.indices = map(lambda x: x + 1, map(int, sMDA.indices))

    def toSelectionMDA(self, reader):

        sMDA = None

        for pr in int2range(self.indices).split():
            if sMDA is None:
                sMDA = reader.parse_selection('bynum ' + pr)
            else:
                sMDA += reader.parse_selection('bynum ' + pr)

        return sMDA
