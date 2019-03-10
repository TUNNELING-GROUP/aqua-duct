# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2018-2019  Tomasz Magdziarz <info@aquaduct.pl>
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

'''
Module performs HDR 2D calculations only with Gaussian Kerneld Density Estimator
as impelemented in :mod:`scipy.stats`.
'''

from aquaduct import logger

iscontour = True
try:
    from matplotlib import _contour
except ImportError:
    logger.warning("Cannot import _contour from matplotlib, contour calculations will be skipped.")
    iscontour = False


def hdr2contour(hdr, fraction=0.9):
    # X,Y,Z,no mask,corner mask,nchunk = 0
    _c = _contour.QuadContourGenerator(hdr.X, hdr.Y, hdr.Z, None, True, 0)
    cc = _c.create_contour(hdr.Z.max() * (1 - fraction))
    if len(cc) == 0:
        return
    return hdr.pca.undo(cc[0], pc=[0, 1])
