# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2017  Tomasz Magdziarz <info@aquaduct.pl>
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
import cPickle as pickle
from importlib import import_module
import json


from aquaduct.utils.helpers import SmartRange,SmartRangeIncrement


################################################################################

class GlobalConfigStore(object):
    cachedir = None
    cachemem = False

GCS = GlobalConfigStore()

class CoordsRangeIndexCache(object):
    cache = {}


CRIC = CoordsRangeIndexCache()

################################################################################
# CIRC save in cache dir

def CRIC_getstate():
    return ((int(k),((int(kk),map(int,SmartRange.raw2sequence(SmartRange(fast_raw=vv.collection).raw_increment))) for kk,vv in v.iteritems())) for k,v in CRIC.cache.iteritems())

def CRIC_setstate(state):
    for k,v in state:
        pass


import collections
import itertools


class FramesRangeCollection(object):
    # currently it is assumed that samrt ranges increments only are possible
    def __init__(self):
        self.collection = [] # order on this list does matter!

    def append(self,srange):
        if not len(self.collection):
            self.collection.append(srange)
            logger.debug("FRC append first srange %s",str(srange))
            return
        # there are V cases:
        #           |------|            sr
        # 1 |---|                       it is before sr
        # 2      |-----|                it overlaps with sr but begins before
        # 3          |----|             it is contained in sr
        # 4              |----|         it overlaps with sr but ends after
        # 24     |------------|         it overlabs with sr in 2 and 4 way
        # 5                  |----|     it is after sr
        while (srange is not None):
            for nr,sr in enumerate(self.collection):
                # sr
                if sr.overlaps_mutual(srange):# or srange.overlaps(sr):
                    if sr.contains(srange):
                        # case 3
                        srange = None
                        break
                    if srange.first_element() < sr.first_element():
                        # case 2
                        self.collection.insert(nr,SmartRangeIncrement(srange.first_element(),sr.first_element()-srange.first_element()))
                        logger.debug("FRC case 2 insert srange %s at %d",str(SmartRangeIncrement(srange.first_element(),sr.first_element()-srange.first_element())),nr)
                        if srange.last_element() > sr.last_element():
                            # case 24
                            srange = SmartRangeIncrement(sr.last_element()+1,srange.last_element()-sr.last_element())
                            break
                        else:
                            srange = None
                        break
                    if srange.last_element() > sr.last_element():
                        # case 4
                        srange = SmartRangeIncrement(sr.last_element()+1,srange.last_element()-sr.last_element())
                        continue
                else:
                    if srange.last_element() < sr.first_element():
                        # case 1: insert it before sr
                        self.collection.insert(nr,srange)
                        logger.debug("FRC case 1 insert srange %s at %d",str(srange),nr)
                        srange = None
                        break
                    if srange.first_element() > sr.last_element():
                        # case 5: do nothing
                        continue
            # if something is left append it to the end
            if srange is not None and nr == len(self.collection) - 1:
                self.collection.append(srange)
                logger.debug("FRC append remaining srage %s",str(srange))
                srange = None

    def get_ranges(self,srange):
        # yield sranges from collection and appropriate ranges for these sranges
        # assumes append was already called? call it!
        # after it is called only case 3 or 4 is possible or no overlap at all
        self.append(srange)
        for sr in self.collection:
            if sr.overlaps(srange):
                if sr.contains(srange):
                    # case 3
                    yield sr,xrange(srange.first_element()-sr.first_element(),srange.first_element()-sr.first_element()+len(srange))
                    srange = None
                    break
                # case 4
                yield sr,xrange(srange.first_element()-sr.first_element(),srange.first_element()-sr.first_element()+sr.last_element()-srange.first_element()+1)
                srange = SmartRangeIncrement(sr.last_element()+1,srange.last_element()-sr.last_element())





#https://stackoverflow.com/questions/12670395/json-encoding-very-long-iterators/45482776#45482776
class IterEncoder(json.JSONEncoder):
    """
    JSON Encoder that encodes iterators as well.
    Write directly to file to use minimal memory
    """
    class FakeListIterator(list):
        def __init__(self, iterable):
            self.iterable = iter(iterable)
            try:
                self.firstitem = next(self.iterable)
                self.truthy = True
            except StopIteration:
                self.truthy = False

        def __iter__(self):
            if not self.truthy:
                return iter([])
            return itertools.chain([self.firstitem], self.iterable)

        def __len__(self):
            raise NotImplementedError("Fakelist has no length")

        def __getitem__(self, i):
            raise NotImplementedError("Fakelist has no getitem")

        def __setitem__(self, i):
            raise NotImplementedError("Fakelist has no setitem")

        def __bool__(self):
            return self.truthy

    def default(self, o):
        if isinstance(o, collections.Iterable):
            return type(self).FakeListIterator(o)
        return super(IterEncoder,self).default(o)


def get_cric_reader(mode='r'):
    # returns file object
    if GCS.cachedir:
        data_file_name = GCS.cachedir + os.path.sep + 'cric.dump' # TODO: magic constant
        if mode=='r' and not os.path.exists(data_file_name):
            return
        if mode=='w' and not os.path.exists(GCS.cachedir):
            os.makedirs(GCS.cachedir)
        logger.debug("Preparing CRIC store with file %s",data_file_name)
        return gzip.open(data_file_name, mode=mode, compresslevel=9)

def save_cric():
    vda = get_cric_reader(mode='w')
    if vda:
        logger.debug("Saving CRIC data.")
        json.dump(CRIC.cache,vda,cls=IterEncoder)
        vda.close()

def load_cric():
    vda = get_cric_reader(mode='r')
    if vda:
        logger.debug("Loading CRIC data.")
        CRIC_state = json.load(vda)
        #CRIC.cache = pickle.load(vda)
        vda.close()
    else:
        CRIC.cache = {}





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

