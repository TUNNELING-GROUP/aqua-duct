# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2017  Tomasz Magdziarz, Alicja Płuciennik, Michał Stolarczyk <info@aquaduct.pl>
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
#

'''
Module implements AutoBarber generation of spheres.
'''


import logging
logger = logging.getLogger(__name__)

from collections import namedtuple

from scipy.spatial.distance import cdist
import numpy as np

from aquaduct.utils import clui
from aquaduct.traj.reader import atom2vdw_radius
from aquaduct.utils.helpers import lind

class Sphere(namedtuple('Sphere', 'center radius')):
    '''
    Simple sphere class.
    '''
    def is_point_within(self,point):
        return self.radius > cdist(np.matrix(self.center), np.matrix(point), metric='euclidean')

    def is_sphere_within(self,sphere):
        center,radius = sphere
        return self.radius > cdist(np.matrix(self.center), np.matrix(center), metric='euclidean') + radius

    def is_sphere_cloud(self,sphere):
        center,radius = sphere
        return self.radius > cdist(np.matrix(self.center), np.matrix(center), metric='euclidean') - radius

class WhereToCut(object):
    '''
    Class implements method for creating (optimal) set of AutoBarber spheres for a collection of spaths;
    access to trajectory is also required to read VdW radii.
    '''
    # creates collection of Spheres
    def __init__(self,spaths,traj_reader,
                 expected_nr_of_spaths=None,
                 selection=None,
                 mincut=None,
                 mincut_level=False,
                 maxcut=None,
                 maxcut_level=False,
                 tovdw=False,
                 forceempty=False):
        '''
        :param list spaths: List of :class:`aquaduct.traj.paths.SinglePath` objects.
        :param traj_reader: :class:`aquaduct.traj.reader.Reader` object.
        :param int expected_nr_of_spaths: Number of spaths passed as :arg:`spaths`. Requilred when length of :arg:`spaths` cannod be calculated, eg when it is a generator..
        :param str selection: Selection string of molecular object used for spheres generation.
        :param float mincut: Value of *mincut* parameter.
        :param float maxcut: Value of *maxcut* parameter.
        :param bool mincut_level: Flag of *mincut_level*.
        :param bool maxcut_level: Flag of *maxcut_level*.
        :param bool tovdw: Flag of to VdW radii correction parameter.
        :param bool forceemtpy: If set *True* spheres of radius 0 are returned if no other sphere can be generated.

        '''
        # force empty means that empty spheres are also returned with radius 0
        self.forceempty = forceempty
        self.selection = selection
        self.mincut = mincut
        self.maxcut = maxcut
        self.mincut_level = mincut_level
        self.maxcut_level = maxcut_level
        self.tovdw = tovdw

        self.expected_nr_of_spaths = expected_nr_of_spaths

        self.spheres = []

        self.add_spheres(spaths, traj_reader)

    def check_minmaxcuts(self):
        mincut = False
        mincut_val = None
        if self.mincut is not None:
            mincut = True
            mincut_val = float(self.mincut)
            logger.debug('mincut set to %0.2f', mincut_val)
        maxcut = False
        maxcut_val = None
        if self.maxcut is not None:
            maxcut = True
            maxcut_val = float(self.maxcut)
            logger.debug('maxcut set to %0.2f', maxcut_val)
        if mincut and maxcut:
            if maxcut_val < mincut_val:
                logger.warning('Values of mincut %0.2f and maxcut %0.2f are mutually exclusive. No spheres will be used in Auto Barber.',mincut_val,maxcut_val)
                if self.maxcut_level or self.mincut_level:
                    logger.warning('Options mincut_level and maxcut_level are set to %s and %s accordingly, some spheres might be retained.')
        return mincut,mincut_val,maxcut,maxcut_val

    def add_spheres(self, spaths, traj_reader):
        clui.message("Auto Barber is looking where to cut:")
        mincut, mincut_val, maxcut, maxcut_val = self.check_minmaxcuts()
        if self.expected_nr_of_spaths:
            pbar = clui.pbar(self.expected_nr_of_spaths)
        else:
            pbar = clui.pbar(len(spaths))
        barber = traj_reader.parse_selection(self.selection)
        vdwradius = 0
        for sp in spaths:
            centers = []
            frames = []
            if sp.has_in:
                centers.append(sp.coords_first_in)
                frames.append(sp.paths_first_in)
            if sp.has_out:
                centers.append(sp.coords_last_out)
                frames.append(sp.paths_last_out)
            for center, frame in zip(centers, frames):
                make_sphere = True
                if make_sphere:
                    traj_reader.set_current_frame(frame)
                    distances = cdist(np.matrix(center), barber.atom_positions(), metric='euclidean').flatten()
                    if self.tovdw:
                        vdwradius = atom2vdw_radius(barber.atoms[np.argmin(distances)])
                    radius = min(distances) - vdwradius
                    if radius <= 0:
                        logger.debug('VdW correction resulted in <= 0 radius.')
                        make_sphere = False
                    if mincut and radius < mincut_val:
                        if not self.mincut_level:
                            logger.debug('Sphere radius %0.2f is less then mincut %0.2f', radius, mincut_val)
                            make_sphere = False
                        else:
                            logger.debug('Sphere radius %0.2f leveled to mincut %0.2f', radius, mincut_val)
                            radius = mincut_val
                    if maxcut and radius > maxcut_val:
                        if not self.maxcut_level:
                            logger.debug('Sphere radius %0.2f is greater then maxcut %0.2f', radius, maxcut_val)
                            make_sphere = False
                        else:
                            logger.debug('Sphere radius %0.2f leveled to maxcut %0.2f', radius, maxcut_val)
                            radius = maxcut_val
                if make_sphere:
                    logger.debug('Added sphere of radius %0.2f' % radius)
                    self.spheres.append(Sphere(center, radius))
                elif self.forceempty:
                    logger.debug('Added sphere of radius 0')
                    self.spheres.append(Sphere(center, 0))
            pbar.next()
        pbar.finish()

    def _cut_thyself(self,spheres_passed):
        if self.expected_nr_of_spaths:
            N = self.expected_nr_of_spaths
        else:
            N = len(spaths)
        # make a deep copy?
        spheres =copy.copy(spheres_passed)
        noredundat_spheres_count = 0
        while True:
            spehres.sort(key=lambda s: s.radius, reverse=True)
            spheres_coords = np.array([sphe.center for sphe in spheres])
            spheres_radii = np.array([sphe.radius for sphe in spheres])
            noredundat_spheres = []
            while spheres:
                big = spheres.pop(0)
                center, radius = big
                distances = cdist(np.matrix(center), spheres_coords[1:], metric='euclidean').flatten()
                # add radii
                distances += spheres_radii[1:]
                # remove if distance <= radius
                to_keep = distances > radius
                # do keep
                spheres_coords = spheres_coords[1:][to_keep]
                spheres_radii = spheres_radii[1:][to_keep]
                # do keep spheres
                to_keep_ids = np.argwhere(to_keep).flatten().tolist()
                spheres = lind(spheres, to_keep_ids)
                # add big to non redundant
                noredundat_spheres.append(big)
            if len(noredundat_spheres) == noredundat_spheres_count:
                break
            else:
                noredundat_spheres_count = len(noredundat_spheres)
                spheres = noredundat_spheres



    def cut_thyself(self):
        # this changes order and removes some of the spheres!
        clui.message("Barber, cut thyself:")
        N = len(self.spheres)
        pbar = clui.pbar(N)
        noredundat_spheres_count = 0
        while True:
            self.spheres.sort(key=lambda s: s.radius, reverse=True)
            # create matrix of coords and vector of radii
            spheres_coords = np.array([sphe.center for sphe in self.spheres])
            spheres_radii = np.array([sphe.radius for sphe in self.spheres])
            noredundat_spheres = []
            while self.spheres:
                # topmost is the biggest one
                big = self.spheres.pop(0)
                center, radius = big
                # calculate distances of all to big but big
                distances = cdist(np.matrix(center), spheres_coords[1:], metric='euclidean').flatten()
                # add radii
                distances += spheres_radii[1:]
                # remove if distance <= radius
                to_keep = distances > radius
                # do keep
                spheres_coords = spheres_coords[1:][to_keep]
                spheres_radii = spheres_radii[1:][to_keep]
                # do keep spheres
                np.argwhere(to_keep).flatten()
                to_keep_ids = np.argwhere(to_keep).flatten().tolist()
                self.spheres = lind(self.spheres, to_keep_ids)
                # add big to non redundant
                noredundat_spheres.append(big)
                pbar.update(N-len(self.spheres))
                logger.debug("Removal of redundant cutting places: done %d, to analyze %d" % (
                    len(noredundat_spheres), len(self.spheres)))
            if len(noredundat_spheres) == noredundat_spheres_count:
                logger.debug("Removal of redundant cutting places done. %d non redundant spheres found." % len(
                    noredundat_spheres))
                break
            else:
                noredundat_spheres_count = len(noredundat_spheres)
                self.spheres = noredundat_spheres
        self.spheres = noredundat_spheres
        pbar.finish()

    def cloud_groups(self):
        # sort reverse
        sortids = np.argsort([sphe.radius for sphe in self.spheres]).tolist()[::-1]
        # start from biggest one
        for spheid in sortids:
            pass


        pass
