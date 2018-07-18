# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2018  Tomasz Magdziarz <info@aquaduct.pl>
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

import cPickle as pickle
import gzip
import os
from collections import OrderedDict

import numpy as np

from aquaduct import version


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

def get_vda_reader(filename, mode="r"):
    if os.path.splitext(filename)[-1].lower() in ['.dump']:
        return ValveDataAccess_pickle(mode=mode, data_file_name=filename)
    elif os.path.splitext(filename)[-1].lower() in ['.npaq']:
        return ValveDataAccess_numpy(mode=mode, data_file_name=filename)
    elif os.path.splitext(filename)[-1].lower() in ['.nc', '.aqnc']:
        return ValveDataAccess_nc(mode=mode, data_file_name=filename)
    raise ValueError('Unknown file type of %s file' % filename)


################################################################################

class ValveDataAccess(object):
    unknown_names = 'UNK'

    def __init__(self, mode=None, data_file_name=None):
        logger.debug('Opening file %s in mode %s' % (data_file_name, mode))
        self.data_file_name = data_file_name
        self.data_file = None
        self.data = None
        assert mode in 'rw'
        self.mode = mode  # r, w
        self.open()

    def open(self):
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

class ValveDataAccess_pickle(ValveDataAccess):
    mimic_old_var_name = 'aq_data_to_save'
    unknown_names = 'UNK'

    def open(self):
        # open file
        self.data_file = gzip.open(self.data_file_name, mode=self.mode, compresslevel=9)
        # if mode is w save header with version etc
        if self.mode == 'w':
            pickle.dump({'version': version(),
                         'aquaduct_version': version()}, self.data_file)
        elif self.mode == 'r':
            versions = pickle.load(self.data_file)
            check_versions(versions)
            # loaded data!
            self.data = {}
            try:
                while True:
                    loaded_data = pickle.load(self.data_file)
                    self.data.update(loaded_data)
            except:  # TODO: remove it!
                pass

    def close(self):
        self.data_file.close()

    def load(self):
        data = {}
        # this is to mimic v0.3 behaviour
        for _name, _value in self.data.iteritems():
            if _name == self.mimic_old_var_name:
                for name, value in _value.iteritems():
                    '''
                    if isinstance(value, CompactSelectionMDA):
                        value = value.toSelectionMDA(self.reader)
                        # with self.reader.get() as traj_reader:
                        #    value = value.toSelectionMDA(traj_reader)
                    '''
                    # TODO: following is to overcome problems with missing names when data is <0.4
                    ################################################################################
                    if name == 'paths':
                        for path_name, path in value.iteritems():
                            if not hasattr(path, 'name'):
                                path.name = self.unknown_names
                    if name == 'spaths':
                        for spath in value:
                            if not hasattr(spath.id, 'name'):
                                spath.id.name = self.unknown_names
                    if name == 'inls':
                        if not hasattr(value, 'spheres'):
                            value.spheres = []
                        if hasattr(value, 'refs'):
                            for r in value.refs:
                                if not hasattr(r, 'name'):
                                    r.name = self.unknown_names
                    ################################################################################
                    data.update({name: value})
                break
            else:
                for name, value in self.data.iteritems():
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
        for name, value in kwargs.iteritems():
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
        if False:  # hasattr(value,'simple_dump'):
            pickle.dump({name: value.simple_dump()}, self.data_file)
        else:
            pickle.dump({name: value}, self.data_file)


ValveDataAccess = ValveDataAccess_pickle
