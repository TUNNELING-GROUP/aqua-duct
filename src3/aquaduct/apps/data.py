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

from aquaduct.utils.helpers import SmartRangeIncrement

logger = logging.getLogger(__name__)

import os
import pickle as pickle
import gzip
import numpy as np
from collections import OrderedDict
from importlib import import_module

#from netCDF4 import Dataset

from aquaduct import version, version_nice, logger
#from aquaduct.traj.selections import CompactSelectionMDA, SelectionMDA
from aquaduct.utils import clui

class GlobalConfigStore(object):
    cachedir = None
    cachemem = False

GCS = GlobalConfigStore()

class CoordsRangeIndexCache(object):
    cache = {}

    def get_frc(self,number,rid):
        # wrapper for get ranges from frc
        if number not in self.cache:
            self.cache.update({number: {}})
            logger.debug("CRIC new number %d", number)
        if rid not in self.cache[number]:
            self.cache[number].update({rid: FramesRangeCollection()})
            logger.debug("CRIC new rid %d", rid)
        return self.cache[number][rid]

    def update_cric(self,cric):
        for number,rid_frc in cric.cache.items():
            for rid,frc in rid_frc.items():
                this_frc = self.get_frc(number,rid)
                for srange in frc.collection:
                    this_frc.append(srange)


CRIC = CoordsRangeIndexCache()

################################################################################
# CIRC save in cache dir

def get_cric_reader(mode='r'):
    if GCS.cachedir:
        try:
            data_file_name = GCS.cachedir + os.path.sep + 'cric.dump'
            if mode=='r' and not os.path.exists(data_file_name):
                return
            if mode=='w' and not os.path.exists(GCS.cachedir):
                os.makedirs(GCS.cachedir)
            vda = get_vda_reader(data_file_name)
            logger.debug("Preparing CRIC store with file %s",data_file_name)
            return vda(mode=mode,data_file_name=data_file_name)
        except IOError:
            logger.warning("Unable to access CRIC data in cache dir [%s]." % GCS.cachedir)
            pass

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
# FRC

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
                    yield sr,range(srange.first_element()-sr.first_element(),srange.first_element()-sr.first_element()+len(srange))
                    srange = None
                    break
                # case 4
                yield sr,range(srange.first_element()-sr.first_element(),srange.first_element()-sr.first_element()+sr.last_element()-srange.first_element()+1)
                srange = SmartRangeIncrement(sr.last_element()+1,srange.last_element()-sr.last_element())


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

'''
class LoadDumpWrapper(object):
    """This is wrapper for pickled data that provides compatibility
    with earlier versions of Aqua-Duct.

    Conversions in use:

    1) replace 'aquaduct.' by 'aquaduct.'

    """

    def __init__(self, filehandle):
        self.fh = filehandle

    def convert(self, s):
        new_s = s
        new_s = new_s.replace('aqueduct.', 'aquaduct.')
        new_s = new_s.replace('aqueduct_version', 'aquaduct_version')
        return new_s

    def read(self, *args, **kwargs):
        return self.convert(self.fh.read(*args, **kwargs))

    def readline(self, *args, **kwargs):
        return self.convert(self.fh.readline(*args, **kwargs))
'''

################################################################################
# VDA reader

def get_vda_reader(filename):
    if os.path.splitext(filename)[-1].lower() in ['.dump']:
        return ValveDataAccess_pickle
    elif os.path.splitext(filename)[-1].lower() in ['.nc','.aqnc']:
        return ValveDataAccess_nc
    raise ValueError('Unknown file type of %s file' % filename)

################################################################################
# Data Access interface

class ValveDataAccess(object):

    def __init__(self, mode=None, data_file_name=None, reader=None):
        logger.debug('Opening file %s in mode %s' % (data_file_name, mode))
        self.data_file_name = data_file_name
        self.data_file = None
        self.data = None
        assert mode in 'rw'
        self.mode = mode  # r, w
        self.reader = reader
        self.open(self.data_file_name, self.mode)

    def open(self, data_file_name, mode):
        raise NotImplementedError()

    def close(self):
        raise NotImplementedError()

    def __del__(self):
        self.close()

    def load(self):
        raise NotImplementedError()

    def dump(self, **kwargs):
        raise NotImplementedError()


################################################################################
# Pickle way

class ValveDataAccess_pickle(ValveDataAccess):
    mimic_old_var_name = 'aq_data_to_save'
    unknown_names = 'UNK'

    def open(self, data_file_name, mode):
        # open file
        self.data_file = gzip.open(data_file_name, mode=mode, compresslevel=9)
        # if mode is w save header with version etc
        if mode == 'w':
            pickle.dump({'version': version(),
                         'aquaduct_version': version()}, self.data_file)
        elif mode == 'r':
            versions = pickle.load(self.data_file)
            check_versions(versions)
            # loaded data!
            self.data = {}
            try:
                while True:
                    loaded_data = pickle.load(self.data_file)
                    self.data.update(loaded_data)
            except:
                pass

    def close(self):
        self.data_file.close()

    def load(self):
        data = {}
        # this is to mimic v0.3 behaviour
        for _name, _value in self.data.items():
            if _name == self.mimic_old_var_name:
                for name, value in _value.items():
                    '''
                    if isinstance(value, CompactSelectionMDA):
                        value = value.toSelectionMDA(self.reader)
                        # with self.reader.get() as traj_reader:
                        #    value = value.toSelectionMDA(traj_reader)
                    '''
                    # TODO: following is to overcome problems with missing names when data is <0.4
                    ################################################################################
                    if name == 'paths':
                        for path_name,path in value.items():
                            if not hasattr(path,'name'):
                                path.name = self.unknown_names
                    if name == 'spaths':
                        for spath in value:
                            if not hasattr(spath.id,'name'):
                                spath.id.name = self.unknown_names
                    if name == 'inls':
                        if not hasattr(value,'spheres'):
                            value.spheres = []
                        if hasattr(value,'refs'):
                            for r in value.refs:
                                if not hasattr(r, 'name'):
                                    r.name = self.unknown_names
                    ################################################################################
                    data.update({name: value})
                break
            else:
                for name, value in self.data.items():
                    '''
                    if isinstance(value, CompactSelectionMDA):
                        value = value.toSelectionMDA(self.reader)
                        # with self.reader.get() as traj_reader:
                        #    value = value.toSelectionMDA(traj_reader)
                    '''
                    # TODO: following is to overcome problems with missing names when data is <0.4
                    ################################################################################
                    '''
                    if name == 'paths':
                        for path_name,path in value.iteritems():
                            if not hasattr(path,'name'):
                                path.name = self.unknown_names
                    if name == 'spaths':
                        for spath in value:
                            if not hasattr(spath.id,'name'):
                                spath.id.name = self.unknown_names
                    if name == 'inls':
                        if not hasattr(value,'spheres'):
                            value.spheres = []
                        if hasattr(value,'refs'):
                            for r in value.refs:
                                if not hasattr(r, 'name'):
                                    r.name = self.unknown_names
                    '''
                    ################################################################################
                    data.update({name: value})
                break
        return data

    def dump(self, **kwargs):
        for name, value in kwargs.items():
            self.set_variable(name, value)

    def get_variable(self, name):
        assert self.mode == "r"
        value = self.data[name]


        return value

    def set_variable(self, name, value):
        assert self.mode == "w"
        '''
        if isinstance(value, SelectionMDA):
            value = CompactSelectionMDA(value)
        '''
        # options are passed as dictionaries already
        '''
        if 'options' in name:
            value = value._asdict()
        '''
        if False: #hasattr(value,'simple_dump'):
            pickle.dump({name: value.simple_dump()}, self.data_file)
        else:
            pickle.dump({name: value}, self.data_file)


################################################################################
# NetCDF helpers
'''
class ValveDataAccessRoots(object):
    roots = []

    def open(self, data_file_name, mode):
        self.roots.append(Dataset(data_file_name, mode, format="NETCDF4"))
        return self.roots[-1]

    def close_all(self):
        for root in self.roots:
            root.close()

    def __del__(self):
        self.close_all()


VDAR = ValveDataAccessRoots()

'''
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
        for k, v in d.items():
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

'''

class ValveDataAccess_nc(ValveDataAccess):

    def __init__(self,*args,**kwargs):
        super(self,ValveDataAccess_nc).__init__(*args,**kwargs)
        self.root = None

    def open(self, data_file_name, mode):
        self.root = VDAR.open(data_file_name, mode)
        logger.debug('Dataset created')

'''

################################################################################
# default way

ValveDataAccess = ValveDataAccess_pickle


