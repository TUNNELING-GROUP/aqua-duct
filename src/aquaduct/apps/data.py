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

import cPickle as pickle
import gzip
import numpy as np
from collections import OrderedDict
from importlib import import_module

from netCDF4 import Dataset

from aquaduct import version, version_nice, logger
from aquaduct.traj.selections import CompactSelectionMDA, SelectionMDA
from aquaduct.utils import clui

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


class ValveDataAccess_pickle(object):
    def __init__(self, mode=None, data_file_name=None, reader=None):
        # as for now. lets use dump files
        self.data_file_name = data_file_name
        self.data_file = None
        self.data = None
        assert mode in 'rw'
        self.mode = mode  # r, w
        self.reader = reader
        self.open(self.data_file_name, self.mode)

    def open(self, data_file_name, mode):
        # open file
        self.data_file = gzip.open(data_file_name, mode=mode, compresslevel=9)
        # if mode is w save header with version etc
        if mode == 'w':
            pickle.dump({'version': version(),
                         'aquaduct_version': version()}, self.data_file)
        elif mode == 'r':
            data_file = LoadDumpWrapper(self.data_file)
            versions = pickle.load(data_file)
            check_versions(versions)
            # loaded data!
            self.data = {}
            try:
                while True:
                    loaded_data = pickle.load(data_file)
                    self.data.update(loaded_data)
            except:
                pass

    def close(self):
        self.data_file.close()

    def __del__(self):
        self.close()

    def get_variable(self, name):
        assert self.mode == "r"
        value = self.data[name]
        if isinstance(value, CompactSelectionMDA):
            with self.reader.get() as traj_reader:
                return value.toSelectionMDA(traj_reader)
        return value

    def load(self):
        data = {}
        for name, value in self.data.iteritems():
            if isinstance(value, CompactSelectionMDA):
                with self.reader.get() as traj_reader:
                    value = value.toSelectionMDA(traj_reader)
            data.update({name: value})
        return data

    def set_variable(self, name, value):
        assert self.mode == "w"
        if isinstance(value, SelectionMDA):
            value = CompactSelectionMDA(value)
        if 'options' in name:
            value = value._asdict()
        pickle.dump({name: value}, self.data_file)

    def dump(self, **kwargs):
        for name, value in kwargs.iteritems():
            self.set_variable(name, value)


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

def get_object_name(something):
    name_ = something.__module__
    if hasattr(something,'__name__'):
        name_ += '.' + something.__name__
    elif hasattr(something,'__class__'):
        if hasattr(something.__class__, '__name__'):
            name_ += '.' + something.__class__.__name__
    return name_

def get_object_from_name(name):
    module_name = '.'.join(name.split('.')[:-1])
    object_name = name.split('.')[-1]
    module = import_module(module_name)
    return getattr(module, object_name)

class ValveDataAccess_nc(object):
    aqt_types_dict_base = OrderedDict(
        {'np': 1, 'list': 2, 'dict': 3, 'dict_int': 32, 'str': 4, 'object': 5, 'bool': 6, 'none': 7,
         'int': 8, 'float': 9,'tuple':16})
    # types names
    aqt_np = 'np'
    aqt_list = 'list'
    aqt_tuple = 'tuple'
    aqt_dict = 'dict'
    aqt_dict_int = 'dict_int'
    aqt_object = 'object'
    aqt_str = 'str'
    aqt_bool = 'bool'
    aqt_none = 'none'
    aqt_int = 'int'
    aqt_float = 'float'

    # metadata names
    aqt_name = 'aqt'
    aqt_object_name = 'aqt_obj_name'
    aqt_object_state = 'aqt_obj_state'
    aqt_types_dict_name = 'aqt_types_dict'
    aqt_dim_name = 'aqt_dim'

    def __init__(self, mode=None, data_file_name=None, reader=None):
        logger.debug('Opening file %s in mode %s' % (data_file_name,mode))
        self.data_file_name = data_file_name
        self.data = None
        assert mode in 'rw'
        self.mode = mode  # r, w
        self.reader = reader
        self.root = VDAR.open(self.data_file_name, self.mode)
        logger.debug('Dataset created')
        self.aqt_types_dict = None
        self.aqt_types_enum = None
        if self.mode == 'w':
            self.aqt_types_enum = self.root.createEnumType(np.int8, self.aqt_types_dict_name, self.aqt_types_dict_base)
        elif self.mode == 'r':
            self.aqt_types_enum = self.root.enumtypes[self.aqt_types_dict_name]
        self.aqt_types_dict = dict(self.aqt_types_enum.enum_dict)

    def create_aqtype(self, gr, aqtype):
        gr.createDimension(self.aqt_dim_name, 1)
        var = gr.createVariable(self.aqt_name, self.aqt_types_enum, (self.aqt_dim_name,))
        var[0] = self.aqt_types_dict[aqtype]
        var = gr.variables[self.aqt_name]
        if var[0] == 0:
            pass

    def get_aqtype(self, gr):
        if self.aqt_name in gr.variables:
            var = gr.variables[self.aqt_name]
            #print var
            #print var[0],self.aqt_types_dict
            if var[0] == 0:
                pass
            return self.aqt_types_dict.keys()[self.aqt_types_dict.values().index(var[0])]
        return self.aqt_dict
        # raise ValueError('Wrong structure')

    def get_dimensions_names(self, name, shape):
        for nr, d in enumerate(shape):
            yield name + str(nr)

    def create_dimensions(self, name, shape, gr):
        out = []
        names = gr.dimensions.keys()
        sizes = [gr.dimensions[dim].size for dim in names]
        for nr, d in enumerate(shape):
            # check if there is such dimension
            if d in sizes:
                out.append(names[sizes.index(d)])
            else:
                gr.createDimension(name + str(nr), d)
                out.append(name + str(nr))
                names.append(name + str(nr))
                sizes.append(d)
        return out

    def create_variable_np(self, name, value, gr):
        names = self.create_dimensions(name, value.shape, gr)
        #print value.dtype
        if value.size > 100:
            var = gr.createVariable(name, value.dtype, tuple(names),zlib=True,shuffle=True,complevel=9)
        else:
            var = gr.createVariable(name, value.dtype, tuple(names))#, zlib=True, shuffle=True, complevel=9)
        var[:] = value

    def set_object(self, name, value, group=None):
        if group is None:
            group = self.root
        # ndarray, create variable
        #print name, type(value), isinstance(value, np.ndarray)
        logger.debug('%s/%s' % (group.path.rstrip('/'), name))
        if isinstance(value, np.ndarray):
            # set variable
            self.create_variable_np(name, value, group)
        elif value is None:
            new_group = group.createGroup(name)
            self.create_aqtype(new_group, self.aqt_none)
            self.set_object('0', np.array(0), group=new_group)  # FIXME: magic constant!
        elif isinstance(value, int):
            new_group = group.createGroup(name)
            self.create_aqtype(new_group, self.aqt_int)
            self.set_object('0', np.array(int(value)), group=new_group)  # FIXME: magic constant!
        elif isinstance(value, float):
            new_group = group.createGroup(name)
            self.create_aqtype(new_group, self.aqt_float)
            self.set_object('0', np.array(float(value)), group=new_group)  # FIXME: magic constant!
        elif isinstance(value, bool):
            new_group = group.createGroup(name)
            self.create_aqtype(new_group, self.aqt_bool)
            self.set_object('0', np.array(int(value)), group=new_group)  # FIXME: magic constant!
        elif isinstance(value, str):
            # make string
            new_group = group.createGroup(name)
            self.create_aqtype(new_group, self.aqt_str)
            self.set_object('0', np.array(value), group=new_group)  # FIXME: magic constant!
        elif isinstance(value, list):
            # create new group
            new_group = group.createGroup(name)
            # make it group of list type
            self.create_aqtype(new_group, self.aqt_list)
            # iterate over the list
            for nr, deep_value in enumerate(value):
                # recurence
                new_name = str(nr)
                self.set_object(new_name, deep_value, new_group)
        elif isinstance(value, tuple):
            # create new group
            new_group = group.createGroup(name)
            # make it group of list type
            self.create_aqtype(new_group, self.aqt_tuple)
            # iterate over the tuple
            for nr, deep_value in enumerate(value):
                # recurence
                new_name = str(nr)
                self.set_object(new_name, deep_value, new_group)
        elif isinstance(value, dict):
            empty = False
            if len(value) == 0:
                empty = True
            # create new group
            new_group = group.createGroup(name)
            # are keys int?
            int_keys = True
            for k in value.keys():
                if not isinstance(k, int):
                    int_keys = False
                    break
            # make it group of list type
            if int_keys:
                self.create_aqtype(new_group, self.aqt_dict_int)
            else:
                self.create_aqtype(new_group, self.aqt_dict)
            # iterate over the list
            for key, deep_value in value.iteritems():
                assert not empty
                # recurence
                new_name = str(key)
                self.set_object(new_name, deep_value, new_group)
        elif isinstance(value, object):
            # create new group
            new_group = group.createGroup(name)
            # make it group of list type
            self.create_aqtype(new_group, self.aqt_object)
            # get full object name
            object_name = get_object_name(value)
            # put marker
            self.set_object(self.aqt_object_name, np.array(object_name), group=new_group)
            # get state
            try:
                object_state = value.__getstate__()
            except AttributeError:
                raise NotImplementedError('Saving is possible for objects implementing __getstate__ method only.')
            # set state
            self.set_object(self.aqt_object_state, object_state, group=new_group)

    def get_object(self, group=None):
        if group is None:
            group = self.root
        aqtype = self.get_aqtype(group)
        logger.debug(group.path)
        if aqtype == self.aqt_dict or aqtype == self.aqt_dict_int:
            out = {}
            # iterate over groups
            for name, new_group in group.groups.iteritems():
                if aqtype == self.aqt_dict_int:
                    out.update({int(name): self.get_object(group=new_group)})
                else:
                    out.update({str(name): self.get_object(group=new_group)})
            for name, var in group.variables.iteritems():
                if name == self.aqt_name: continue
                if aqtype == self.aqt_dict_int:
                    out.update({int(name): var[:]})
                else:
                    out.update({str(name): var[:]})
        elif aqtype in [self.aqt_list, self.aqt_tuple]:
            ids = group.groups.keys()
            ids += group.variables.keys()
            ids.sort()
            out = []
            for idd in ids:
                if idd == self.aqt_name: continue
                if idd in group.groups:
                    out.append(self.get_object(group=group.groups[idd]))
                elif idd in group.variables:
                    out.append(group.variables[idd][:])
            if aqtype == self.aqt_tuple:
                out = tuple(out)
        elif aqtype == self.aqt_str:
            var = group.variables['0']
            out = str(var[:])
        elif aqtype == self.aqt_bool:
            var = group.variables['0']
            out = bool(var[:])
        elif aqtype == self.aqt_int:
            var = group.variables['0']
            out = int(var[:])
        elif aqtype == self.aqt_float:
            var = group.variables['0']
            out = float(var[:])
        elif aqtype == self.aqt_none:
            out = None
        elif aqtype == self.aqt_object:
            object_name = group.variables[self.aqt_object_name][:]
            state = self.get_object(group=group.groups[self.aqt_object_state])
            #print object_name, state
            object_instance = get_object_from_name(object_name)
            object_instance = object_instance.__new__(object_instance)
            #print type(object_instance)
            object_instance.__setstate__(state, reader=self.reader)
            out = object_instance

        return out

    def dump(self, **kwargs):
        self.set_object(name='root', value=kwargs)

    def load(self):
        return self.get_object()['root']


# default
ValveDataAccess = ValveDataAccess_nc
