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
        pickle.dump(CRIC.cache,vda)
        vda.close()

def load_cric():
    vda = get_cric_reader(mode='r')
    if vda:
        logger.debug("Loading CRIC data.")
        CRIC.cache = pickle.load(vda)
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

