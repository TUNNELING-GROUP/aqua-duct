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

import logging

logger = logging.getLogger(__name__)

from itertools import izip_longest

import numpy as np

import MDAnalysis as mda

from aquaduct.geom.convexhull import SciPyConvexHull
from aquaduct.geom.convexhull import is_point_within_convexhull
from aquaduct.utils.helpers import int2range, are_rows_uniq
from aquaduct.utils.maths import make_default_array

class Selection(object):
    """
    def __init__(self,selection,selection_string=None):

        self.selection_object = selection
        self.selection_string = selection_string
    """

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
        return SciPyConvexHull(self.atom_positions())

    def contains_residues(self,other,convex_hull=False,map_fun=None):
        # Checks if this selection contains other.
        # If convex_hull is False it checks if residues of this selection matches to other.
        # If convex_hull is True it checks if residues of other selection are within convex hull of this.
        # Returns list of bool values.
        # Parameter map_fun is optional map function. If it is None, regular map will be used.
        if map_fun is None:
            map_fun = map
        if convex_hull:
            other_coords = list(other.center_of_mass_of_residues())
            chull = self.get_convexhull_of_atom_positions()
            return map_fun(is_point_within_convexhull, izip_longest(other_coords, [], fillvalue=chull))
        else:
            # check if other selection is empty
            if other.unique_resids_number() == 0:
                return []
            this_uids = self.unique_resids(ikwid=True)
            return [res_other.unique_resids(ikwid=True) in this_uids for res_other in other.iterate_over_residues()]

    def containing_residues(self,other,*args,**kwargs):
        # Convienience wrapper for contains_residues.
        #def get_res_in_scope(is_res_in_scope, res):
        other_new = None
        for iris, other_residue in zip(self.contains_residues(other, *args, **kwargs), other.iterate_over_residues()):
            if iris:
                if other_new is None:
                    other_new = other_residue
                else:
                    other_new += other_residue
        return other_new


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
        # check if positions are correct
        if len(self.atoms):
            positions = self.positions
            if not are_rows_uniq(positions):
                logger.warning('Some atoms have the same positions, your data might me corrupted.')
            return make_default_array(positions)
        logger.warning('Selection comprises no atoms, check your settings.')
        return make_default_array([])

    def __add__(self, other):
        if other is not None:
            return SelectionMDA(self._atoms + other._atoms)
        return SelectionMDA(self._atoms)

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
