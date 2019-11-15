# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2017-2019  Tomasz Magdziarz <info@aquaduct.pl>
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

import os
import gzip
from importlib import import_module
import json
import collections

from aquaduct.utils.helpers import SmartRange, SmartRangeIncrement, SmartRangeEqual


################################################################################

class GlobalConfigStore(object):
    _cachedir = None
    _cachemem = False
    netcdf = False
    nc4 = False
    sandwich_import = False
    cachetype = 'full'

    @property
    def cachedir(self):
        return self._cachedir

    @cachedir.setter
    def cachedir(self, value):
        assert not self.sandwich_import, "Cachedir property cannot be set after sandiwch import."
        self._cachedir = value

    @property
    def cachemem(self):
        return self._cachemem

    @cachemem.setter
    def cachemem(self, value):
        assert not self.sandwich_import, "Cachemem property cannot be set after sandiwch import."
        self._cachemem = value


GCS = GlobalConfigStore()


class CoordsRangeIndexCache(object):
    cache = {}

    def get_frc(self, number, rid):
        # wrapper for get ranges from frc
        if number not in self.cache:
            self.cache.update({number: {}})
            logger.debug("CRIC new number %d", number)
        if rid not in self.cache[number]:
            self.cache[number].update({rid: FramesRangeCollection()})
            logger.debug("CRIC new rid %d", rid)
        return self.cache[number][rid]

    def update_cric(self, cric):
        for number, rid_frc in cric.cache.items():
            for rid, frc in rid_frc.items():
                this_frc = self.get_frc(number, rid)
                for srange in frc.collection:
                    this_frc.append(srange)

    def setstate(self, state):

        for number, rid_frc in state:
            for rid, frc in rid_frc:
                this_frc = self.get_frc(number, rid)
                for srange in SmartRange(fast_minc_seq=frc).raw:
                    this_frc.append(srange)

    def reset(self):
        self.cache = {}

    def getstate(self):
        # this is saved as json so it constitutes format
        return ((int(k),
                 ((int(kk), list(map(int, SmartRange.raw2sequence(SmartRange(fast_raw=vv.collection).raw_increment)))) for
                  kk, vv
                  in v.items())) for k, v in self.cache.items())


CRIC = CoordsRangeIndexCache()


################################################################################
# CIRC save in cache dir

def get_cric_reader(mode='r'):
    # returns file object
    if GCS.cachedir:
        try:
            data_file_name = GCS.cachedir + os.path.sep + 'cric.json'  # TODO: magic constant
            if mode == 'r' and not os.path.exists(data_file_name):
                return
            if mode == 'w' and not os.path.exists(GCS.cachedir):
                os.makedirs(GCS.cachedir)
            logger.debug("Preparing CRIC store with file %s", data_file_name)
            return gzip.open(data_file_name, mode=mode+'t', compresslevel=6)
        except IOError:
            logger.warning("Unable to access CRIC data in cache dir [%s]." % GCS.cachedir)
            pass


def save_cric():
    vda = get_cric_reader(mode='w')
    if vda:
        logger.debug("Saving CRIC data.")
        json.dump(CRIC.getstate(), vda, cls=IterEncoder, indent=None, separators=(',', ':'))
        vda.close()


def load_cric():
    vda = get_cric_reader(mode='r')
    if vda:
        logger.debug("Loading CRIC data.")
        CRIC.reset()
        CRIC.setstate(json.load(vda))
        vda.close()
    else:
        CRIC.reset()


################################################################################
# FRC

class FramesRangeCollection(object):
    # currently it is assumed that samrt ranges increments only are possible
    def __init__(self):
        self.collection = []  # order on this list does matter!

    def append(self, srange):
        # TODO: remove it later
        if isinstance(srange, SmartRangeEqual):
            srange = SmartRangeIncrement(srange.element, srange.times)
        if not len(self.collection):
            self.collection.append(srange)
            logger.debug("FRC append first srange %s", str(srange))
            return
        # there are V cases:
        #           |------|            sr
        # 1 |---|                       it is before sr
        # 2      |-----|                it overlaps with sr but begins before
        # 3          |----|             it is contained in sr
        # 4              |----|         it overlaps with sr but ends after
        # 24     |------------|         it overlabs with sr in 2 and 4 way
        # 5                  |----|     it is after sr
        while srange is not None:
            for nr, sr in enumerate(self.collection):
                # sr
                if sr.overlaps_mutual(srange):  # or srange.overlaps(sr):
                    if sr.contains(srange):
                        # case 3
                        srange = None
                        break
                    if srange.first_element() < sr.first_element():
                        # case 2
                        self.collection.insert(nr, SmartRangeIncrement(srange.first_element(),
                                                                       sr.first_element() - srange.first_element()))
                        logger.debug("FRC case 2 insert srange %s at %d", str(
                            SmartRangeIncrement(srange.first_element(), sr.first_element() - srange.first_element())),
                                     nr)
                        if srange.last_element() > sr.last_element():
                            # case 24
                            srange = SmartRangeIncrement(sr.last_element() + 1,
                                                         srange.last_element() - sr.last_element())
                            break
                        else:
                            srange = None
                        break
                    if srange.last_element() > sr.last_element():
                        # case 4
                        srange = SmartRangeIncrement(sr.last_element() + 1, srange.last_element() - sr.last_element())
                        continue
                else:
                    if srange.last_element() < sr.first_element():
                        # case 1: insert it before sr
                        self.collection.insert(nr, srange)
                        logger.debug("FRC case 1 insert srange %s at %d", str(srange), nr)
                        srange = None
                        break
                    if srange.first_element() > sr.last_element():
                        # case 5: do nothing
                        continue
            # if something is left append it to the end
            if srange is not None and nr == len(self.collection) - 1:
                self.collection.append(srange)
                logger.debug("FRC append remaining srange %s", str(srange))
                srange = None

    def get_ranges(self, srange):
        # yield sranges from collection and appropriate ranges for these sranges
        # assumes append was already called? call it!
        # after it is called only case 3 or 4 is possible or no overlap at all
        self.append(srange)
        for sr in self.collection:
            if sr.overlaps(srange):
                # TODO: remove it later
                if isinstance(sr, SmartRangeEqual):
                    sr = SmartRangeIncrement(sr.element, sr.times)
                if sr.contains(srange):
                    # case 3
                    yield sr, range(srange.first_element() - sr.first_element(),
                                     srange.first_element() - sr.first_element() + len(srange))
                    srange = None
                    break
                # case 4
                yield sr, range(srange.first_element() - sr.first_element(),
                                 srange.first_element() - sr.first_element() +
                                 sr.last_element() - srange.first_element() + 1)
                srange = SmartRangeIncrement(sr.last_element() + 1,
                                             srange.last_element() - sr.last_element())


class IterEncoder(json.JSONEncoder):

    def default(self, o):
        if isinstance(o, collections.Iterable):
            return list(o)
        return super(IterEncoder, self).default(o)


################################################################################


def get_object_name(something):
    name_ = something.__module__
    if hasattr(something, '__name__'):
        name_ += '.' + something.__name__
    elif hasattr(something, '__class__'):
        if hasattr(something.__class__, '__name__'):
            name_ += '.' + something.__class__.__name__
    return name_


def get_object_from_name(name):
    module_name = '.'.join(name.split('.')[:-1])
    object_name = name.split('.')[-1]
    module = import_module(module_name)
    return getattr(module, object_name)

################################################################################
