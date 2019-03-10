# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2016-2018  Tomasz Magdziarz, Alicja Płuciennik, Michał Stolarczyk <info@aquaduct.pl>
# Copyright (C) 2019  Tomasz Magdziarz <info@aquaduct.pl>
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


from aquaduct import logger

from collections import namedtuple

import numpy as np
from scipy.spatial.distance import cdist
import copy

from aquaduct.utils import clui
from aquaduct.utils.helpers import lind
from aquaduct import __mail__


class Sphere(namedtuple('Sphere', 'center radius nr')):
    """
    Simple sphere class.
    """

    def is_point_within(self, point):
        return self.radius > cdist(np.matrix(self.center), np.matrix(point), metric='euclidean')

    def is_sphere_within(self, sphere):
        center, radius = sphere
        return self.radius > cdist(np.matrix(self.center), np.matrix(center), metric='euclidean') + radius

    def is_sphere_cloud(self, sphere):
        center, radius = sphere
        return self.radius > cdist(np.matrix(self.center), np.matrix(center), metric='euclidean') - radius


def do_cut_thyself(spheres_passed, progress=False):
    # returns noredundant spheres
    # make a deep copy?
    # TODO: this is not memory efficient?
    spheres = copy.copy(spheres_passed)
    N = len(spheres)
    if progress:
        clui.message("Barber, cut thyself:")
        pbar = clui.pbar(N)
    noredundat_spheres_count = 0
    redundat_spheres = []
    while True:
        spheres.sort(key=lambda s: s.radius, reverse=True)
        spheres_coords = np.array([sphe.center for sphe in spheres])
        spheres_radii = np.array([sphe.radius for sphe in spheres])
        noredundat_spheres = []
        while spheres:
            big = spheres.pop(0)
            center, radius, nr = big
            distances = cdist(np.matrix(center), spheres_coords[1:], metric='euclidean').flatten()
            # add radii
            distances += spheres_radii[1:]
            # remove if distance <= radius
            to_keep = distances > radius
            to_remove = ~to_keep
            # do keep
            spheres_coords = spheres_coords[1:][to_keep]
            spheres_radii = spheres_radii[1:][to_keep]
            # do keep spheres
            to_keep_ids = np.argwhere(to_keep).flatten().tolist()
            to_remove_ids = np.argwhere(to_remove).flatten().tolist()
            redundat_spheres.extend(lind(spheres, to_remove_ids))
            spheres = lind(spheres, to_keep_ids)
            # add big to non redundant
            noredundat_spheres.append(big)
            if progress:
                # pbar.update(N - len(spheres))
                pbar.update(len(redundat_spheres))
            logger.debug("Removal of redundant cutting places: done %d, to analyze %d" % (
                len(noredundat_spheres), len(spheres)))
        if len(noredundat_spheres) == noredundat_spheres_count:
            logger.debug("Removal of redundant cutting places done. %d non redundant spheres found." % len(
                noredundat_spheres))
            break
        else:
            noredundat_spheres_count = len(noredundat_spheres)
            spheres = noredundat_spheres
    if progress:
        pbar.finish()

    assert len(noredundat_spheres) + len(
        redundat_spheres) == N, "Inconsistent number of not and redundant spheres. Please send a bug report to the developer(s): %s" % __mail__
    # sorting
    noredundat_spheres.sort(key=lambda s: s.nr)
    redundat_spheres.sort(key=lambda s: s.nr)
    return noredundat_spheres, redundat_spheres
