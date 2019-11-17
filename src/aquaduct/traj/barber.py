# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2017-2018  Tomasz Magdziarz, Alicja Płuciennik, Michał Stolarczyk <info@aquaduct.pl>
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
#

"""
Module implements AutoBarber generation of spheres.
"""

from aquaduct import logger
from scipy.spatial.distance import cdist
import numpy as np

from aquaduct.utils import clui
from aquaduct.utils.helpers import listify
from aquaduct.utils.helpers import lind, SmartRange
from aquaduct.traj.sandwich import ReaderAccess, Reader



from aquaduct.geom import Sphere, do_cut_thyself
from aquaduct.utils.multip import optimal_threads
from multiprocessing import Pool
from functools import partial
from itertools import chain
from aquaduct.apps.data import CRIC

__mail__ = 'info@aquaduct.pl'


@listify
def spaths2spheres(spaths, minmax=None, selection=None, tovdw=None, forceempty=None):
    mincut, mincut_val, maxcut, maxcut_val, mincut_level, maxcut_level = minmax

    for sp in spaths:

        traj_reader = Reader.get_reader_by_id(sp.id.id).open()
        barber = traj_reader.parse_selection(selection)

        vdwradius = 0

        centers = []
        frames = []
        # TODO: This is inconsistent with inlets types. Below is equivalent to surface only inlets. Rework this and make it coherent to each other.
        # Assume it is already coherent.
        if sp.has_in:
            centers.append(sp.coords_first_in)
            frames.append(sp.paths_first_in)
        if sp.has_out:
            centers.append(sp.coords_last_out)
            frames.append(sp.paths_last_out)
        for center, frame in zip(centers, frames):
            make_sphere = True
            if make_sphere:
                traj_reader.set_frame(frame)
                distances = cdist(np.matrix(center), np.matrix(list(barber.coords())), metric='euclidean').flatten()
                if tovdw:
                    vdwradius = list(barber.ix(np.argmin(distances)).vdw())[0]
                    logger.debug('VdW correction %0.2f', vdwradius)
                radius = min(distances) - vdwradius
                if radius <= 0:
                    logger.debug('VdW correction resulted in <= 0 radius.')
                    make_sphere = False
                if mincut and radius < mincut_val:
                    if not mincut_level:
                        logger.debug('Sphere radius %0.2f is less then mincut %0.2f', radius, mincut_val)
                        make_sphere = False
                    else:
                        logger.debug('Sphere radius %0.2f leveled to mincut %0.2f', radius, mincut_val)
                        radius = mincut_val
                if maxcut and radius > maxcut_val:
                    if not maxcut_level:
                        logger.debug('Sphere radius %0.2f is greater then maxcut %0.2f', radius, maxcut_val)
                        make_sphere = False
                    else:
                        logger.debug('Sphere radius %0.2f leveled to maxcut %0.2f', radius, maxcut_val)
                        radius = maxcut_val
            if make_sphere:
                logger.debug('Added sphere of radius %0.2f' % radius)
                yield (center.flatten(), radius)
            elif forceempty:
                logger.debug('Added sphere of radius 0')
                yield (center.flatten(), 0)
    yield CRIC


@listify
def inlets2spheres(inlets, minmax=None, selection=None, tovdw=None, forceempty=None):
    mincut, mincut_val, maxcut, maxcut_val, mincut_level, maxcut_level = minmax

    for inlet in inlets:

        traj_reader = Reader.get_reader_by_id(inlet.reference.id).open()
        barber = traj_reader.parse_selection(selection)

        vdwradius = 0

        center = inlet.coords
        frame = inlet.frame

        make_sphere = True
        if make_sphere:
            traj_reader.set_frame(frame)
            distances = cdist(np.matrix(center), np.matrix(list(barber.coords())), metric='euclidean').flatten()
            if tovdw:
                vdwradius = list(barber.ix(np.argmin(distances)).vdw())[0]
                logger.debug('VdW correction %0.2f', vdwradius)
            radius = min(distances) - vdwradius
            if radius <= 0:
                logger.debug('VdW correction resulted in <= 0 radius.')
                make_sphere = False
            if mincut and radius < mincut_val:
                if not mincut_level:
                    logger.debug('Sphere radius %0.2f is less then mincut %0.2f', radius, mincut_val)
                    make_sphere = False
                else:
                    logger.debug('Sphere radius %0.2f leveled to mincut %0.2f', radius, mincut_val)
                    radius = mincut_val
            if maxcut and radius > maxcut_val:
                if not maxcut_level:
                    logger.debug('Sphere radius %0.2f is greater then maxcut %0.2f', radius, maxcut_val)
                    make_sphere = False
                else:
                    logger.debug('Sphere radius %0.2f leveled to maxcut %0.2f', radius, maxcut_val)
                    radius = maxcut_val
        if make_sphere:
            logger.debug('Added sphere of radius %0.2f' % radius)
            yield (center.flatten(), radius)
        elif forceempty:
            logger.debug('Added sphere of radius 0')
            yield (center.flatten(), 0)
    yield CRIC


class WhereToCut(ReaderAccess):
    '''
    Class implements method for creating (optimal) set of AutoBarber spheres for a collection of spaths;
    access to trajectory is also required to read VdW radii.
    '''

    # creates collection of Spheres
    def __init__(self,
                 spaths=None,
                 inlets=None,
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
        :param int expected_nr_of_spaths: Number of spaths passed. Requilred when length of spaths cannod be calculated, eg when it is a generator.
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

        if spaths is not None:
            self.add_spheres_from_spaths(spaths)
        if inlets is not None:
            self.add_spheres_from_inlets(inlets)

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
                logger.warning(
                    'Values of mincut %0.2f and maxcut %0.2f are mutually exclusive. No spheres will be used in Auto Barber.',
                    mincut_val, maxcut_val)
                if self.maxcut_level or self.mincut_level:
                    logger.warning(
                        'Options mincut_level and maxcut_level are set to %s and %s accordingly, some spheres might be retained.')
        return mincut, mincut_val, maxcut, maxcut_val

    def add_spheres_from_spaths(self, spaths):
        clui.message("Auto Barber is looking where to cut:")
        if self.expected_nr_of_spaths:
            pbar = clui.pbar(self.expected_nr_of_spaths)
            n_add = self.expected_nr_of_spaths
        else:
            pbar = clui.pbar(len(spaths))
            n_add = len(spaths)
        minmax = self.check_minmaxcuts() + (self.mincut_level, self.maxcut_level)


        Reader.reset()
        pool = Pool(processes=optimal_threads.threads_count)
        n = max(1, optimal_threads.threads_count)
        chunks = (n if chunk <= (n_add / n - 1) else (n_add % n) for chunk in range(int(n_add / n) + int(np.sign(n_add % n))))

        add_function = partial(spaths2spheres, minmax=minmax, selection=self.selection, tovdw=self.tovdw,
                               forceempty=self.forceempty)
        _spaths = chain(spaths)
        Reader.reset()
        spheres_new = pool.imap(add_function, ([next(_spaths) for cc in range(c)] for c in chunks))

        nr = 0
        for spheres in spheres_new:
            for center_radius_or_cric in spheres:
                if isinstance(center_radius_or_cric, tuple):
                    center, radius = center_radius_or_cric
                    self.spheres.append(Sphere(center=center, radius=radius, nr=self.get_current_nr()))
                else:
                    CRIC.update_cric(center_radius_or_cric)
            nr += n
            if nr < n_add:
                pbar.update(nr)
            else:
                pbar.update(n_add)
        pool.close()
        pool.join()
        pbar.finish()

    def add_spheres_from_inlets(self, inlets):
        clui.message("Auto Barber is looking where to cut:")
        if self.expected_nr_of_spaths:
            pbar = clui.pbar(self.expected_nr_of_spaths)
            n_add = self.expected_nr_of_spaths
        else:
            pbar = clui.pbar(len(inlets.inlets_list))
            n_add = len(inlets.inlets_list)
        minmax = self.check_minmaxcuts() + (self.mincut_level, self.maxcut_level)


        Reader.reset()
        pool = Pool(processes=optimal_threads.threads_count)
        # pool = Pool(processes=1)
        n = max(1, optimal_threads.threads_count)
        chunks = (n if chunk <= (n_add / n - 1) else (n_add % n) for chunk in range(n_add / n + np.sign(n_add % n)))

        add_function = partial(inlets2spheres, minmax=minmax, selection=self.selection, tovdw=self.tovdw,
                               forceempty=self.forceempty)
        _inlets = chain(inlets)
        spheres_new = pool.imap(add_function, [[next(_inlets) for cc in range(c)] for c in chunks])
        # spheres_new = imap(add_function, ([_inlets.next() for cc in xrange(c)] for c in chunks))

        nr = 0
        for spheres in spheres_new:
            for center_radius_or_cric in spheres:
                if isinstance(center_radius_or_cric, tuple):
                    center, radius = center_radius_or_cric
                    self.spheres.append(Sphere(center=center, radius=radius, nr=self.get_current_nr()))
                else:
                    CRIC.update_cric(center_radius_or_cric)
            nr += n
            if nr < n_add:
                pbar.update(nr)
            else:
                pbar.update(n_add)
        pool.close()
        pool.join()
        pbar.finish()

    def get_current_nr(self):
        if len(self.spheres):
            # nr = max((sphe.nr for sphe in self.spheres))
            nr = self.spheres[-1].nr
            nr += 1
        else:
            nr = 0
        assert nr >= len(self.spheres), "Inconsistent number of spheres."
        return nr

    def inlet2sphere(self, inlet):


        traj_reader = Reader.get_reader_by_id(inlet.reference.id).open()
        mincut, mincut_val, maxcut, maxcut_val = self.check_minmaxcuts()
        barber = traj_reader.parse_selection(self.selection)
        vdwradius = 0

        center = inlet.coords
        frame = inlet.frame

        make_sphere = True
        if make_sphere:
            traj_reader.set_frame(frame)
            distances = cdist(np.matrix(center), np.matrix(list(barber.coords())), metric='euclidean').flatten()
            if self.tovdw:
                vdwradius = list(barber.ix(np.argmin(distances)).vdw())[0]
                logger.debug('VdW correction %0.2f', vdwradius)
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
            return Sphere(center.flatten(), radius, self.get_current_nr())
        elif self.forceempty:
            logger.debug('Added sphere of radius 0')
            return Sphere(center.flatten(), 0, self.get_current_nr())

    def spath2spheres(self, sp):

        traj_reader = self.get_reader_by_id(sp.id.id)
        mincut, mincut_val, maxcut, maxcut_val = self.check_minmaxcuts()
        barber = traj_reader.parse_selection(self.selection)
        vdwradius = 0

        centers = []
        frames = []
        # TODO: This is inconsistent with inlets types. Below is equivalent to surface only inlets. Rework this and make it coherent to each other.
        # Assume it is already coherent.
        if sp.has_in:
            centers.append(sp.coords_first_in)
            frames.append(sp.paths_first_in)
        if sp.has_out:
            centers.append(sp.coords_last_out)
            frames.append(sp.paths_last_out)
        for center, frame in zip(centers, frames):
            make_sphere = True
            if make_sphere:
                traj_reader.set_frame(frame)
                distances = cdist(np.matrix(center), np.matrix(list(barber.coords())), metric='euclidean').flatten()
                if self.tovdw:
                    vdwradius = list(barber.ix(np.argmin(distances)).vdw())[0]
                    logger.debug('VdW correction %0.2f', vdwradius)
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
                yield Sphere(center.flatten(), radius, self.get_current_nr())
            elif self.forceempty:
                logger.debug('Added sphere of radius 0')
                yield Sphere(center.flatten(), 0, self.get_current_nr())

    def cut_thyself(self):
        self.spheres = do_cut_thyself(self.spheres, progress=True)[0]

    def is_overlaping_with_cloud(self, sphere):
        spheres_coords = np.array([sphe.center for sphe in self.spheres])
        spheres_radii = np.array([sphe.radius for sphe in self.spheres])
        center, radius, nr = sphere
        distances = cdist(np.matrix(center), spheres_coords, metric='euclidean').flatten()
        distances = distances - spheres_radii - radius
        return (distances <= 0).any()

    def cloud_groups(self, progress=False):
        # no redundant spheres
        noredundant_spheres, redundant_spheres = do_cut_thyself(self.spheres, progress=progress)
        if progress:
            clui.message("Barber, clouds clusters:")
            pbar = clui.pbar(len(self.spheres))

        noredundant_spheres_coords = np.array([sphe.center for sphe in noredundant_spheres])
        noredundant_spheres_radii = np.array([sphe.radius for sphe in noredundant_spheres])
        clouds = {}
        for nrs in noredundant_spheres:
            # calculate distance
            center, radius, nr = nrs
            distances = cdist(np.matrix(center), noredundant_spheres_coords, metric='euclidean').flatten()
            distances = distances - noredundant_spheres_radii - radius
            current_cloud = set(np.argwhere(distances <= 0).flatten().tolist())
            del distances
            # check if cci overlaps with any of already found clouds
            cloud_id_intersections = []
            for cloud_id, cloud in clouds.items():
                if current_cloud.intersection(cloud):
                    # current cloud intersects with cloud
                    cloud_id_intersections.append(cloud_id)
            # does it intersects?
            if cloud_id_intersections:
                for cii in cloud_id_intersections:
                    current_cloud = current_cloud.union(clouds.pop(cii))
                    # current id?
            current_id = list(clouds.keys())
            if current_id:
                for cid in range(max(current_id) + 2):
                    if cid not in current_id:
                        current_id = cid
                        break
            else:
                current_id = 0
            clouds.update({current_id: current_cloud})
            if progress:
                pbar.next()

        # chnage nrs id to global ids; add redundant spheres
        nrs_gids = [nrs.nr for nrs in noredundant_spheres]
        for cloud_id, cloud in clouds.items():

            cloud = sorted(list(cloud))

            if len(redundant_spheres):
                cloud_coords = noredundant_spheres_coords[cloud]
                cloud_radii = noredundant_spheres_radii[cloud]

            cloud = sorted([nrs_gids[c] for c in cloud])
            # clouds.update({cloud_id:cloud})

            for rs_id in range(len(redundant_spheres))[::-1]:
                center, radius, nr = redundant_spheres[rs_id]
                distances = cdist(np.matrix(center), cloud_coords, metric='euclidean').flatten()
                distances = distances - cloud_radii - radius
                current_cloud = np.argwhere(distances <= 0).flatten().tolist()
                if current_cloud:
                    cloud.append(nr)
                    redundant_spheres.pop(rs_id)
                    if progress:
                        pbar.next()

            clouds.update({cloud_id: cloud})

        if progress:
            pbar.finish()

        return clouds


def barber_with_spheres(coords, spheres):
    '''
    param numpy.ndarray coords: Path's coordinates subjected to barber procedure.
    param spheres: Spheres used to cut input coordinates.
    rtype: :class:`numpy.ndarray`.
    return: List of indices to be kept.
    '''
    tokeep = np.ones(len(coords), dtype=np.bool)
    for sp in spheres:
        tokeep = np.logical_and(tokeep, cdist(coords, np.array([sp.center]), metric='euclidean').flatten() > sp.radius)
    return np.argwhere(tokeep).flatten().tolist()


def barber_with_spheres_big_matrix(coords, spheres):
    # calculate big distance matrix
    # distances = cdist(self.coords, [s.center for s in spheres],metric='euclidean')
    # compare with radii
    # tokeep = distances > np.matrix([[s.radius for s in spheres]])
    if len(spheres):
        return np.argwhere((cdist(np.array(list(coords)), [s.center for s in spheres], metric='euclidean') > np.matrix(
            [[s.radius for s in spheres]])).all(1).A1).flatten().tolist()


@listify
def barber_paths(paths, spheres=None, only_for_names=None):
    # cut paths with barber
    for path in (pat for pat in paths if pat.name in only_for_names):
        tokeep = barber_with_spheres(path.coords, spheres)
        path.update_types_frames(SmartRange(lind(path.types, tokeep)), SmartRange(lind(path.frames, tokeep)))
        yield path
    yield CRIC
