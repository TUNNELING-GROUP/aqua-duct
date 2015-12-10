


import MDAnalysis as mda
import numpy as np

#from aquarium.geom.convexhull_pyhull import ConvexHull
from aquarium.geom.convexhull import ConvexHull

class Selection(object):

    '''
    def __init__(self,selection,selection_string=None):

        self.selection_object = selection
        self.selection_string = selection_string
    '''

    def center_of_mass(self):
        # should return numpy (3,) array
        raise NotImplementedError()

    def iterate_over_residues(self):
        # should iterate over residues in the selection returning object of the same type
        raise NotImplementedError()

    def unique_resids(self):
        # should return array of resids
        raise NotImplementedError()

    def atom_positions(self):
        # should return numpy (x,3) array
        raise NotImplementedError()

    def center_of_mass_of_residues(self):
        # should resturn list of lists or generator of center of masses
        return (R.center_of_mass().tolist() for R in self.iterate_over_residues())

    def get_convexhull_of_atom_positions(self):
        # should return modified ConvexHull object
        return ConvexHull(self.atom_positions())

    def uniquify(self):
        # should change selection to unique atoms only
        raise NotImplementedError()

    def __add__(self, other):
        raise NotImplementedError()


class SelectionMDA(mda.core.AtomGroup.AtomGroup,
                   Selection):

    def iterate_over_residues(self):
        return (self.__class__(R) for R in self.residues)

    def unique_resids(self):
        return np.unique(self.resids)

    def atom_positions(self):
        return self.positions

    def __add__(self, other):
        return SelectionMDA(self._atoms + other._atoms)

    def uniquify(self):
        self.__init__(mda.core.AtomGroup.AtomGroup(set(self._atoms)))