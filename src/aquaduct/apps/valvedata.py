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

from scipy.io import netcdf
import numpy as np

from aquaduct import version
from aquaduct.utils.helpers import dictify

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

def get_vda_reader(filename,mode="r"):
    if os.path.splitext(filename)[-1].lower() in ['.dump']:
        return ValveDataAccess_pickle(mode=mode,data_file_name=filename)
    elif os.path.splitext(filename)[-1].lower() in ['.npaq']:
        return ValveDataAccess_numpy(mode=mode, data_file_name=filename)
    elif os.path.splitext(filename)[-1].lower() in ['.nc','.aqnc']:
        return ValveDataAccess_nc(mode=mode,data_file_name=filename)
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

    def get_variable(self, name):
        raise NotImplementedError()

    def set_variable(self, name, value):
        raise NotImplementedError()

################################################################################

from aquaduct.traj.sandwich import ResidueSelection
from itertools import chain
import netCDF4 as nc4


class ValveDataCodec(object):
    # this is in fact definition of data format
    # it assumes data is dictionary with scipy netcdf variables, would it work with netcdf4?

    @staticmethod
    def varname(name,*suffix):
        suff = '.'.join(map(str,suffix))
        return '%s.%s' % (name,suff)

    @staticmethod
    def encode(name,value):
        if name == 'center_of_system':
            yield name,np.array(value)
        if name == 'all_res':
            # save layers
            layers = np.array(value.selected.keys(),dtype=np.int32)
            yield ValveDataCodec.varname(name,'layers'),layers
            # save each layer as separate array
            for l in layers:
                yield ValveDataCodec.varname(name,'layer',l),np.array(value.selected[l],dtype=np.int32)
        if name == 'number_frame_rid_in_object':
            # number of layers
            layers = len(value)
            yield ValveDataCodec.varname(name,'layers'),np.array(layers,dtype=np.int32)
            for nr,layer in enumerate(value):
                # for each layer make sizes array
                yield ValveDataCodec.varname(name,'layer',nr,'sizes'),np.fromiter((len(row) for row in layer),dtype=np.int32)
                # then save layer
                yield ValveDataCodec.varname(name, 'layer', nr), np.fromiter(chain(*(row for row in layer)),dtype=np.int32)
            # make one long array
            # number of layers,
            # size of each layer
            # sizes of all rows
            # rows
            '''
            data_iter = chain([len(value)],
                              (len(layer) for layer in value),
                              chain(*((len(row) for row in layer) for layer in value)),
                              chain(*(chain(*(row for row in layer)) for layer in value)))
            yield name,np.fromiter(data_iter,dtype=np.int32)
            '''

    @staticmethod
    def decode(name,data):
        if name == 'center_of_system':
            return data[name][:].copy()
        if name == 'all_res':
            # create empty all_res object
            # read layers
            layers = data[ValveDataCodec.varname(name,'layers')][:].copy()
            return ResidueSelection(((l,data[ValveDataCodec.varname(name,'layer',l)]) for l in layers))
        if name == 'number_frame_rid_in_object':
            out = []
            # number of layers
            layers = int(data[ValveDataCodec.varname(name,'layers')][:].copy())
            for nr in xrange(layers):
                out_ = []
                sizes = data[ValveDataCodec.varname(name,'layer',nr,'sizes')]
                layer = data[ValveDataCodec.varname(name,'layer',nr)]
                seek = 0
                for s in sizes:
                    out_.append(layer[seek:seek+s].copy())
                    seek += s
                out.append(out_)

            '''                    
            seek = 0
            layers_nr = int(data[name][seek].copy())
            seek += 1
            layers_size = data[name][seek:seek + layers_nr].copy()
            seek += layers_nr
            rows_size = data[name][seek:seek + sum(layers_size)].copy()
            seek += sum(layers_size)
            out = []
            ls_cum = 0
            for l,ls in zip(xrange(layers_nr),layers_size):
                out.append([])
                for rs in rows_size[ls_cum:ls_cum+ls]:
                    row = data[name][seek:seek+rs].copy()
                    seek += rs
                    out[-1].append(list(row))
            '''
            return out



class ValveDataAccess_nc(ValveDataAccess):


    def open(self):
        #self.data_file = netcdf.netcdf_file(self.data_file_name,self.mode)
        self.data_file = nc4.Dataset(self.data_file_name,self.mode)
        if self.mode == 'w':
            self.set_variable('version',np.array(version(),dtype=np.int16))
            self.set_variable('aquaduct_version',np.array(version(),dtype=np.int16))
        elif self.mode == 'r':
            versions =  dict([(k,tuple(v)) for k,v in zip(['version','aquaduct_version'],
                                                          (self.get_variable('version'),
                                                           self.get_variable('aquaduct_version')))])
            check_versions(versions)

    def close(self):
        self.data_file.close()

    def get_variable(self, name,copy=True):
        print name
        if copy:
            return self.data_file.variables[name][:].copy()
        return self.data_file.variables[name]

    def set_variable(self, name, value):
        print name,value
        assert self.mode == "w"
        # value has to be ndarray
        assert isinstance(value,np.ndarray)
        # create dimensions
        dimensions = []
        for nr,d in enumerate(value.shape):
            dimensions.append("%s%d" % (name,nr))
            self.data_file.createDimension(dimensions[-1],d)
        # create variable
        v = self.data_file.createVariable(name,value.dtype,tuple(dimensions),zlib=True)
        # fill variable
        v[:] = value

    def dump(self, **kwargs):
        for name, value in kwargs.iteritems():
            for nname,vvalue in ValveDataCodec.encode(name,value):
                self.set_variable(nname,vvalue)

    @dictify
    def load(self):
        # names of data objects
        names = list(set([name.split('.')[0] for name in self.data_file.variables.keys() if name not in ['version','aquaduct_version']]))
        all_names = [name for name in self.data_file.variables.keys() if name not in ['version','aquaduct_version']]
        for name in names:
            # read all parts of object and decode
            this_object = dict(((n,self.get_variable(n,copy=False)) for n in all_names if name in n))
            yield name,ValveDataCodec.decode(name,this_object)
            for k in this_object.keys():
                v = this_object.pop(k)[:].copy()
                del v

################################################################################

class ValveDataAccess_numpy(ValveDataAccess):
    pass


################################################################################

class ValveDataAccess_pickle(ValveDataAccess):
    mimic_old_var_name = 'aq_data_to_save'

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
            except:
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
        if False: #hasattr(value,'simple_dump'):
            pickle.dump({name: value.simple_dump()}, self.data_file)
        else:
            pickle.dump({name: value}, self.data_file)

