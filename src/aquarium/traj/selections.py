


import MDAnalysis as mda

class Selection(object):

    '''
    def __init__(self,selection,selection_string=None):

        self.selection_object = selection
        self.selection_string = selection_string
    '''

    def center_of_mass(self):
        # should return numpy (3,) array
        raise NotImplementedError()


class SelectionMDA(mda.core.AtomGroup.AtomGroup,
                   Selection):

    pass