


import MDAnalysis as mda
import numpy as np

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

    def unique_resnums(self):
        # should return array of resids
        raise NotImplementedError()


class SelectionMDA(mda.core.AtomGroup.AtomGroup,
                   Selection):


    def iterate_over_residues(self):
        return (self.__class__(R) for R in self.residues)

    def unique_resnums(self):
        return np.unique(self.resids)
