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

from aquaduct.apps.valvedata import ValveDataAccess_pickle, ValveDataAccessRoots, get_vda_reader

logger = logging.getLogger(__name__)

import os
import numpy as np
from collections import OrderedDict
from importlib import import_module

from aquaduct import version, logger
#from aquaduct.traj.selections import CompactSelectionMDA, SelectionMDA

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
    if GCS.cachedir:
        data_file_name = GCS.cachedir + os.path.sep + 'cric.dump'
        if mode=='r' and not os.path.exists(data_file_name):
            return
        if mode=='w' and not os.path.exists(GCS.cachedir):
            os.makedirs(GCS.cachedir)
        vda = get_vda_reader(data_file_name)
        logger.debug("Preparing CRIC store with file %s",data_file_name)
        return vda(mode=mode,data_file_name=data_file_name)

def save_cric():
    vda = get_cric_reader(mode='w')
    if vda:
        logger.debug("Saving CRIC data.")
        vda.dump(**{"CRIC":CRIC.cache})

def load_cric():
    vda = get_cric_reader(mode='r')
    if vda:
        logger.debug("Loading CRIC data.")
        CRIC.cache = vda.load()['CRIC']
    else:
        CRIC.cache = {}


################################################################################
# Version checking

def check_version_compliance(current, loaded, what):
    if current[0] > loaded[0]:
        logger.error('Loaded data has %s major version lower then the application.' % what)
    if current[0] < loaded[0]:
        logger.error('Loaded data has %s major version higher then the application.' % what)
    if current[0] != loaded[0]:
        logger.error('Possible problems with API compliance.')
    if current[1] > loaded[1]:
        logger.warning('Loaded data has %s minor version lower then the application.' % what)
    if current[1] < loaded[1]:
        logger.warning('Loaded data has %s minor version higher then the application.' % what)
    if current[1] != loaded[1]:
        logger.warning('Possible problems with API compliance.')


def check_versions(version_dict):
    assert isinstance(version_dict, (dict, OrderedDict)), "File is corrupted, cannot read version data."
    assert 'version' in version_dict, "File is corrupted, cannot read version data."
    assert 'aquaduct_version' in version_dict, "File is corrupted, cannot read version data."
    check_version_compliance(version(), version_dict['aquaduct_version'], 'Aqua-Duct')
    check_version_compliance(version(), version_dict['version'], 'Valve')


################################################################################
# Pickle compatibility with earlier versions

################################################################################
# VDA reader

################################################################################
# Data Access interface


################################################################################
# numpy way


################################################################################
# Pickle way


################################################################################
# NetCDF helpers


VDAR = ValveDataAccessRoots()


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
# 2array2 objects' decoders

class IdsOverIds(object):
    # use it for netcdf
    @staticmethod
    def dict2arrays(d):
        values = []
        keys_lens = []
        for k, v in d.iteritems():
            keys_lens.append((k, len(v)))
            values.extend(v.tolist())
        return {'values': np.array(values), 'keys_lens': np.array(keys_lens)}

    @staticmethod
    def arrays2dict(values=None, keys_lens=None):
        out = {}
        ll = 0
        for kl in keys_lens:
            k = kl[0]
            l = kl[1]
            v = values[ll:ll + l]
            ll += l
            out.update({int(k): v.tolist()})
        return out


################################################################################
# NetCDF way


################################################################################
# default way

ValveDataAccess = ValveDataAccess_pickle
