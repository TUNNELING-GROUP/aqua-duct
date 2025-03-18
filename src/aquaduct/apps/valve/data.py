# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2018-2019  Tomasz Magdziarz, Michał Banas <info@aquaduct.pl>
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

from aquaduct import logger
# import logging
# logger = logging.getLogger(__name__)

from aquaduct.utils.clui import SimpleTree
import pickle as pickle
import gzip
import os
from collections import OrderedDict

from aquaduct.apps.data import GCS

try:
    if GCS.nc4:
        import netCDF4 as netcdf

        logger.debug('NetCDF AQ format: Imported netCDF4')
    else:
        raise ImportError()
except ImportError:
    from scipy.io import netcdf

    logger.debug('NetCDF AQ format: Imported scipy.io.netcdf')

import numpy as np

from aquaduct import version
from aquaduct.utils.helpers import dictify

from aquaduct.traj.sandwich import ResidueSelection
from aquaduct.traj.paths import GenericPaths, SinglePath, PassingPath, SinglePathID, GenericPathTypeCodes, MasterPath
from aquaduct.utils.helpers import SmartRange
from aquaduct.traj.inlets import InletClusterExtendedType, InletClusterGenericType, Inlet, Inlets
from aquaduct.traj.barber import Sphere
from aquaduct.geom.master import FakeSingleResidueSelection

from itertools import chain
import json


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

    def get_variable(self, name):
        raise NotImplementedError()

    def set_variable(self, name, value):
        raise NotImplementedError()


################################################################################


class ValveDataCodec(object):
    '''
    Class defines format of encoding AQ objects into NetCDF format.
    Encoding defined here can be also used to store AQ objects as NumPy arrays.
    '''

    # this is in fact definition of data format
    # it assumes data is dictionary with scipy netcdf variables, would it work with netcdf4?

    version = 0, 0, 1
    '''
    Current version of :class:`ValveDataCodec`
    '''

    @staticmethod
    def varname(name, *suffix):
        '''
        Name of variable made by combining base name and suffixes (if any).
        Base name and suffixes are joined with dot '.'.

        :param name: Base for variable name.
        :param suffix: Optional suffixes.
        :return: Name of variable made by combining base name and suffixes (if any).
        Base name and suffixes are joined with dot '.'.
        '''
        if len(suffix):
            suff = '.'.join(map(str, suffix))
            return '%s.%s' % (name, suff)
        return '%s' % name

    @staticmethod
    def encode(name, value):

        if name == 'center_of_system':
            # center_of_system: (3,)*float

            yield name, np.array(value)

        if name == 'all_res':
            # all_res.layers: (L,)*int
            #   L is a number of layers
            # all_res.layer.N: (S,)*int
            #   N is a consecutive layer number
            #   S is a number of traced molecules in N layer

            # save layers
            layers = np.array(list(value.selected.keys()), dtype=np.int32)
            yield ValveDataCodec.varname(name, 'layers'), layers
            # save each layer as separate array
            for l in layers:
                yield ValveDataCodec.varname(name, 'layer', l), np.array(value.selected[l], dtype=np.int32)

        if name == 'number_frame_rid_in_object':
            # number_frame_rid_in_object.layers: (1,)*int
            #   number of layers
            # number_frame_rid_in_object.layer.N.sizes: (F,)*int
            #   N is a consecutive layer number
            #   F is a number of frames in layer N
            # number_frame_rid_in_object.layer.N: (Q,)*int
            #   Q is a number of molecules in in object in all frames,
            #   it is a sum of number_frame_rid_in_object.layer.N.sizes

            # TODO It can be further improved because there is no sense to store repeated information for frames.
            # number of layers
            layers = len(value)
            yield ValveDataCodec.varname(name, 'layers'), np.array([layers], dtype=np.int32)
            for nr, layer in enumerate(value):
                # for each layer make sizes array
                yield ValveDataCodec.varname(name, 'layer', nr, 'sizes'), np.fromiter((len(row) for row in layer),
                                                                                      dtype=np.int32)
                # then save layer
                yield ValveDataCodec.varname(name, 'layer', nr), np.fromiter(chain(*(row for row in layer)),
                                                                             dtype=np.int32)

        if name == 'paths':
            # paths.layers:  (L,)*int
            #   L is a number of layers
            # paths.layer.N.names: (P,3)*str
            #   N is a consecutive layer number
            #   P is a number of paths in N layer
            # paths.layer.N.ids: (P,)*int
            #   N is a consecutive layer number
            #   P is a number of paths in N layer
            # paths.layer.N.min_max_frames: (1,2)*int
            # paths.layer.N.object.sizes: (P,)*int
            #   P is a number of paths in N layer
            # paths.layer.N.scope.sizes: (P,)*int
            #   P is a number of paths in N layer
            # paths.layer.N.object: (PO*2,)*int
            #   PO is a number of pairs (start,times) decoding object paths for all paths in N layer
            # paths.layer.N.scope: (PS*2,)*int
            #   PO is a number of pairs (start,times) decoding scope paths for all paths in N layer

            # get all layers
            layers = sorted(set((p.id[0] for p in value)))
            yield ValveDataCodec.varname(name, 'layers'), np.array(layers, dtype=np.int32)
            for N in layers:
                P = sum((1 for p in value if p.id[0] == N))
                yield ValveDataCodec.varname(name, 'layer', N, 'names'), np.fromiter(
                    chain(*(p.name for p in value if p.id[0] == N)), dtype='S1').reshape(P, 3)
                yield ValveDataCodec.varname(name, 'layer', N, 'ids'), np.fromiter(
                    (p.id[-1] for p in value if p.id[0] == N), dtype=np.int32)
                mmf = next((p for p in value if p.id[0] == N))
                yield ValveDataCodec.varname(name, 'layer', N, 'min_max_frames'), np.array(
                    (mmf.min_possible_frame, mmf.max_possible_frame), dtype=np.int32)
                oors = [SmartRange(p.frames_of_object) for p in value if p.id[0] == N]
                yield ValveDataCodec.varname(name, 'layer', N, 'object', 'sizes'), np.fromiter(
                    (2 * len(list(o.raw_increment)) for o in oors), dtype=np.int32)
                yield ValveDataCodec.varname(name, 'layer', N, 'object'), np.fromiter(
                    chain(*(o.raw2sequence(o.raw_increment) for o in oors)), dtype=np.int32)
                oors = [SmartRange(p.frames_of_scope) for p in value if p.id[0] == N]
                yield ValveDataCodec.varname(name, 'layer', N, 'scope', 'sizes'), np.fromiter(
                    (2 * len(list(o.raw_increment)) for o in oors), dtype=np.int32)
                yield ValveDataCodec.varname(name, 'layer', N, 'scope'), np.fromiter(
                    chain(*(o.raw2sequence(o.raw_increment) for o in oors)), dtype=np.int32)

        if name == 'spaths':
            # spaths.layers:  (L,)*int
            #   L is a number of layers
            # spaths.layer.N.names: (P,3)*str
            #   N is a consecutive layer number
            #   P is a number of paths in N layer
            # spaths.layer.N.ids: (P,2)*int
            #   N is a consecutive layer number
            #   P is a number of paths in N layer
            # spaths.layer.N.single: (P,)*int
            #   P is a number of paths in N layer
            # spaths.layer.N.frames: (P,5)*int
            #   P is a number of paths in N layer
            #   table of (begining,end,in,object,out)
            # spaths.layer.N.object.sizes: (PSO,)*int
            #   PSO is a number of paths in N layer that are not PassingPaths
            # spaths.layer.N.object: (PSOS*2,)*int
            #   PSOS is a number of pairs (start,times) decoding strict obejct frames of non PassingPaths for all paths in N layer

            layers = sorted(set((p.id.id[0] for p in value)))
            yield ValveDataCodec.varname(name, 'layers'), np.array(layers, dtype=np.int32)
            for N in layers:
                P = sum((1 for p in value if p.id.id[0] == N))  # number of paths in N
                yield ValveDataCodec.varname(name, 'layer', N, 'names'), np.fromiter(
                    chain(*(p.id.name for p in value if p.id.id[0] == N)), dtype='S1').reshape(P, 3)
                yield ValveDataCodec.varname(name, 'layer', N, 'ids'), np.fromiter(
                    chain(*((p.id.id[-1], p.id.nr) for p in value if p.id.id[0] == N)), dtype=np.int32).reshape(P, 2)
                yield ValveDataCodec.varname(name, 'layer', N, 'single'), np.fromiter(
                    (p.is_single() for p in value if p.id.id[0] == N), dtype='i1')  # i1 instead of bool
                yield ValveDataCodec.varname(name, 'layer', N, 'frames'), np.fromiter(
                    chain(*((p.begins, p.ends) + tuple(p.sizes) for p in value if p.id.id[0] == N)),
                    dtype=np.int32).reshape(P, 5)
                osf = [SmartRange(p.path_object_strict()) for p in value if p.id.id[0] == N and p.is_single()]
                yield ValveDataCodec.varname(name, 'layer', N, 'object', 'sizes'), np.fromiter(
                    (2 * len(list(o.raw_increment)) for o in osf), dtype=np.int32)
                yield ValveDataCodec.varname(name, 'layer', N, 'object'), np.fromiter(
                    chain(*(o.raw2sequence(o.raw_increment) for o in osf)), dtype=np.int32)

        if name == 'ctypes':
            # ctypes: (P,4)*int
            #   None is replaced by -1

            yield ValveDataCodec.varname(name), np.fromiter(
                chain(*([-1 if c is None else c for c in ct.clusters] for ct in value)), dtype=np.int32).reshape(
                len(value), 4)

        if name in ['master_paths', 'master_paths_smooth']:
            # master_paths*.keys: (M,2)*int
            #   M is a number of master paths
            #   keys are ctypes generic
            # master_paths*.names: (M,3)*str
            #   M is a number of master paths
            # master_paths*.ids: (M,3)*int
            #   M is a number of master paths
            # master_paths*.frames: (M,5)*int
            #   M is a number of master paths
            #   table of (begining,end,in,object,out)
            # master_paths*.widths: (MS,)*float
            #   MS is a sum of lenghts of all master paths
            # master_paths*.coords: (MS,3)*int
            #   MS is a sum of lenghts of all master paths
            # master_paths*.object.sizes: (M,)*int
            #   M is a number of master paths
            # master_paths*.object: (MOS*2,)*int
            #   MOS is a number of pairs (start,times) decoding strict obejct frames of master paths

            M = len(value)
            yield ValveDataCodec.varname(name, 'keys'), np.fromiter(
                chain(*([-1 if c is None else c for c in ct.clusters] for ct in list(value.keys()))),
                dtype=np.int32).reshape(len(value), 2)
            yield ValveDataCodec.varname(name, 'names'), np.fromiter(chain(*(p.id.name for p in list(value.values()))),
                                                                     dtype='S1').reshape(M, 3)
            yield ValveDataCodec.varname(name, 'ids'), np.fromiter(
                chain(*((p.id.id[0], p.id.id[-1], p.id.nr) for p in list(value.values()))), dtype=np.int32).reshape(M, 3)
            yield ValveDataCodec.varname(name, 'frames'), np.fromiter(
                chain(*((p.begins, p.ends) + tuple(p.sizes) for p in list(value.values()))), dtype=np.int32).reshape(M, 5)
            yield ValveDataCodec.varname(name, 'widths'), np.fromiter(chain(*(p.width_cont for p in list(value.values()))),
                                                                      dtype=np.float32)
            MS = sum((p.size for p in list(value.values())))
            yield ValveDataCodec.varname(name, 'coords'), np.fromiter(
                chain(*chain(*(p.coords_cont for p in list(value.values())))), dtype=np.float32).reshape(MS, 3)
            osf = [SmartRange(p.path_object_strict()) for p in list(value.values())]
            yield ValveDataCodec.varname(name, 'object', 'sizes'), np.fromiter(
                (2 * len(list(o.raw_increment)) for o in osf), dtype=np.int32)
            yield ValveDataCodec.varname(name, 'object'), np.fromiter(
                chain(*(o.raw2sequence(o.raw_increment) for o in osf)), dtype=np.int32)

        if name == 'inls':
            # inls is a complicated object
            #
            # inls.center_of_system: (3,)*float
            # inls.onlytype: (1,)*str
            # inls.inlets_list.coords: (I,3)*float
            # inls.inlets_list.frame: (I,)*int
            # inls.inlets_list.type: (I,)*int
            #   it references onlytype
            # inls.inlets_list.reference.ids: (I,3)*int
            # inls.inlets_list.reference.name: (I,3)*str
            # inls.inlets_ids: (I,)*int
            #   I is a number of inlets
            # inls.clusters: (I,)*int
            #   I is a number of inlets
            # inls.number_of_clustered_inlets: (1,)*int
            # inls.spheres: (I,4)*float
            #   I is a number of inlets
            # inls.spheres.nr: (I,)*int
            #   I is a number of inlets
            # inls.passing: (1,)*int
            # inls.tree: (1,)*str

            if value.center_of_system is not None:
                yield ValveDataCodec.varname(name, 'center_of_system'), np.array(value.center_of_system)
            # yield ValveDataCodec.varname(name, 'onlytype'), np.fromiter(json.dumps(value.onlytype),dtype='S1')
            yield ValveDataCodec.varname(name, 'inlets_list', 'coords'), np.array([i.coords for i in value.inlets_list])
            yield ValveDataCodec.varname(name, 'inlets_list', 'frame'), np.array([i.frame for i in value.inlets_list],
                                                                                 dtype=np.int32)
            yield ValveDataCodec.varname(name, 'inlets_list', 'type'), np.array(
                [value.onlytype.index(i.type) for i in value.inlets_list], dtype=np.int32)
            yield ValveDataCodec.varname(name, 'inlets_list', 'reference', 'ids'), np.array(
                [list(i.reference.id) + [i.reference.nr] for i in value.inlets_list], dtype=np.int32)
            I = len(value.inlets_ids)
            yield ValveDataCodec.varname(name, 'inlets_list', 'reference', 'name'), np.fromiter(
                chain(*[i.reference.name for i in value.inlets_list]), dtype='S1').reshape(I, 3)
            yield ValveDataCodec.varname(name, 'inlets_ids'), np.array(value.inlets_ids, dtype=np.int32)
            yield ValveDataCodec.varname(name, 'clusters'), np.array(value.clusters, dtype=np.int32)
            yield ValveDataCodec.varname(name, 'number_of_clustered_inlets'), np.array(
                [value.number_of_clustered_inlets],
                dtype=np.int32)
            yield ValveDataCodec.varname(name, 'spheres'), np.array(
                [s.center.tolist() + [s.radius] for s in value.spheres])
            yield ValveDataCodec.varname(name, 'spheres', 'nr'), np.array([s.nr for s in value.spheres], dtype=np.int32)
            yield ValveDataCodec.varname(name, 'passing'), np.array([value.passing], dtype='i1')
            yield ValveDataCodec.varname(name, 'tree'), np.fromiter(repr(value.tree), dtype='S1')

    @staticmethod
    def decode(name, data):
        if name == 'center_of_system':
            return data[name][:].copy()
        if name == 'all_res':
            # create empty all_res object
            # read layers
            layers = data[ValveDataCodec.varname(name, 'layers')][:].copy()
            return ResidueSelection(((l, data[ValveDataCodec.varname(name, 'layer', l)]) for l in layers))
        if name == 'number_frame_rid_in_object':
            out = []
            # number of layers
            layers = int(data[ValveDataCodec.varname(name, 'layers')][:].copy())
            for nr in range(layers):
                out_ = []
                sizes = data[ValveDataCodec.varname(name, 'layer', nr, 'sizes')]
                layer = data[ValveDataCodec.varname(name, 'layer', nr)]
                seek = 0
                for s in sizes:
                    out_.append(layer[seek:seek + s].copy())
                    seek += s
                out.append(out_)
            return out
        if name == 'paths':
            out = []
            layers = data[ValveDataCodec.varname(name, 'layers')]
            for N in layers:
                names = data[ValveDataCodec.varname(name, 'layer', N, 'names')]
                ids = data[ValveDataCodec.varname(name, 'layer', N, 'ids')]
                mmf = data[ValveDataCodec.varname(name, 'layer', N, 'min_max_frames')]
                object_sizes = data[ValveDataCodec.varname(name, 'layer', N, 'object', 'sizes')]
                object_frames = data[ValveDataCodec.varname(name, 'layer', N, 'object')]
                scope_sizes = data[ValveDataCodec.varname(name, 'layer', N, 'scope', 'sizes')]
                scope_frames = data[ValveDataCodec.varname(name, 'layer', N, 'scope')]
                seek_object = 0
                seek_scope = 0
                for osize, ssize, n, pid in zip(object_sizes, scope_sizes, names, ids):
                    out.append(
                        GenericPaths(tuple(map(int, (N, pid))), name_of_res=str(n.tostring()), min_pf=int(mmf[0]),
                                     max_pf=int(mmf[1])))
                    foo = list(SmartRange(fast_minc_seq=object_frames[seek_object:seek_object + osize]).get())
                    seek_object += osize
                    fos = list(SmartRange(fast_minc_seq=scope_frames[seek_scope:seek_scope + ssize]).get())
                    seek_scope += ssize
                    out[-1].add_foos(foo, fos)
            return out
        if name == 'spaths':
            out = []
            layers = data[ValveDataCodec.varname(name, 'layers')]
            for N in layers:
                names = data[ValveDataCodec.varname(name, 'layer', N, 'names')]
                ids = data[ValveDataCodec.varname(name, 'layer', N, 'ids')]
                is_single = data[ValveDataCodec.varname(name, 'layer', N, 'single')]
                frames_table = data[ValveDataCodec.varname(name, 'layer', N, 'frames')]
                object_strict = data[ValveDataCodec.varname(name, 'layer', N, 'object')]
                single_path_nr = 0
                seek_object = 0
                for n, pid, ft, iss in zip(names, ids, frames_table, is_single):
                    spid = SinglePathID(path_id=tuple(map(int, (N, pid[0]))), nr=int(pid[-1]), name=str(n.tostring()))
                    path = list(range(ft[0], ft[1] + 1))
                    if iss:
                        out.append(SinglePath(spid, [[], [], []], [[], [], []]))
                        seek = 0
                        path_in = path[seek:seek + ft[2]]
                        seek += ft[2]
                        path_object = path[seek:seek + ft[3]]
                        seek += ft[3]
                        path_out = path[seek:seek + ft[4]]
                        seek += ft[4]
                        path_object_strict_size = data[ValveDataCodec.varname(name, 'layer', N, 'object', 'sizes')][
                            single_path_nr]
                        single_path_nr += 1
                        path_object_strict = list(SmartRange(
                            fast_minc_seq=object_strict[seek_object:seek_object + path_object_strict_size]).get())
                        seek_object += path_object_strict_size
                        out[-1].add_paths4(path_in, path_object, path_object_strict, path_out)
                    else:
                        out.append(PassingPath(spid, path, [GenericPathTypeCodes.scope_name] * len(path)))
            return out

        if name == 'ctypes':
            ctypes = data[ValveDataCodec.varname(name)]
            return [InletClusterExtendedType(cc[0], cc[2], cc[3], cc[1]) for cc in ([None if c == -1 else c for c in ct] for ct in ctypes)]

        if name in ['master_paths', 'master_paths_smooth']:
            out = {}
            keys = data[ValveDataCodec.varname(name, 'keys')]
            names = data[ValveDataCodec.varname(name, 'names')]
            ids = data[ValveDataCodec.varname(name, 'ids')]
            frames_table = data[ValveDataCodec.varname(name, 'frames')]
            widths = data[ValveDataCodec.varname(name, 'widths')]
            coords = data[ValveDataCodec.varname(name, 'coords')]  # corrds use the same seek as widths
            object_strict = data[ValveDataCodec.varname(name, 'object')]
            seek_widths = 0
            seek_object = 0
            for mpnr, (k, n, pid, ft) in enumerate(zip(keys, names, ids, frames_table)):
                key = InletClusterGenericType(*[None if c == -1 else c for c in k])
                spid = SinglePathID(path_id=tuple(map(int, pid[:2])), nr=int(pid[-1]), name=str(n.tostring()))
                path = list(range(ft[0], ft[1] + 1))

                mp = SinglePath(spid, [[], [], []], [[], [], []])
                seek = 0
                path_in = path[seek:seek + ft[2]]
                seek += ft[2]
                path_object = path[seek:seek + ft[3]]
                seek += ft[3]
                path_out = path[seek:seek + ft[4]]
                seek += ft[4]

                path_object_strict_size = data[ValveDataCodec.varname(name, 'object', 'sizes')][mpnr]
                path_object_strict = list(
                    SmartRange(fast_minc_seq=object_strict[seek_object:seek_object + path_object_strict_size]).get())
                seek_object += path_object_strict_size

                fsrs = FakeSingleResidueSelection(spid.id, range(ft[0], ft[1] + 1),
                                                  coords[seek_widths:seek_widths + ft[1] - ft[0] + 1].copy())
                mp = SinglePath(spid, [[], [], []], [[], [], []], single_res_selection=fsrs)

                mp.add_paths4(path_in, path_object, path_object_strict, path_out)
                mp = MasterPath(mp)
                mp.add_width(widths[seek_widths:seek_widths + ft[1] - ft[0] + 1].copy())
                seek_widths += (ft[1] - ft[0] + 1)

                out.update({key: mp})
            return out

        if name == 'inls':
            if ValveDataCodec.varname(name, 'center_of_system') in data:
                center_of_system = data[ValveDataCodec.varname(name, 'center_of_system')][:].copy()
            else:
                center_of_system = None
            # onlytype = data[ValveDataCodec.varname(name, 'onlytype')][:].copy()
            # onlytype = json.loads(str(onlytype.tostring()))
            # onlytype = [tuple(ot) for ot in onlytype]
            passing = bool(data[ValveDataCodec.varname(name, 'passing')][:].copy())
            inls = Inlets([], center_of_system=center_of_system, passing=passing, onlytype=onlytype)

            coords = data[ValveDataCodec.varname(name, 'inlets_list', 'coords')]
            frame = data[ValveDataCodec.varname(name, 'inlets_list', 'frame')]
            type_ = data[ValveDataCodec.varname(name, 'inlets_list', 'type')]
            ref_ids = data[ValveDataCodec.varname(name, 'inlets_list', 'reference', 'ids')]
            ref_name = data[ValveDataCodec.varname(name, 'inlets_list', 'reference', 'name')]
            for c, f, t, pid, n in zip(coords, frame, type_, ref_ids, ref_name):
                spid = SinglePathID(path_id=tuple(map(int, pid[:2].copy())), nr=int(pid[-1].copy()),
                                    name=str(n.copy().tostring()))
                # inls.inlets_list.append(Inlet(coords=c.copy(), type=onlytype[t.copy()], reference=spid, frame=f.copy()))
                inls.inlets_list.append(Inlet(coords=c.copy(), reference=spid, frame=f.copy()))

            inls.inlets_ids = data[ValveDataCodec.varname(name, 'inlets_ids')][:].copy().tolist()
            inls.clusters = data[ValveDataCodec.varname(name, 'clusters')][:].copy().tolist()
            inls.number_of_clustered_inlets = int(
                data[ValveDataCodec.varname(name, 'number_of_clustered_inlets')][:].copy())

            spheres = data[ValveDataCodec.varname(name, 'spheres')]
            spheres_nr = data[ValveDataCodec.varname(name, 'spheres', 'nr')]
            for s, nr in zip(spheres, spheres_nr):
                inls.spheres.append(Sphere(s[:3].copy(), s[-1].copy(), nr.copy()))

            inls.tree = SimpleTree(treestr=str(data[ValveDataCodec.varname(name, 'tree')][:].copy().tostring()))

            return inls


class ValveDataAccess_nc(ValveDataAccess):
    not_variable = ['version', 'aquaduct_version', 'ValveDataCodec']

    def open(self):
        if hasattr(netcdf, 'Dataset'):
            # netcdf4
            self.data_file = netcdf.Dataset(self.data_file_name, self.mode)
        else:
            # scipy
            self.data_file = netcdf.netcdf_file(self.data_file_name, self.mode)
        if self.mode == 'w':
            self.set_variable('version', np.array(version(), dtype=np.int16))
            self.set_variable('aquaduct_version', np.array(version(), dtype=np.int16))
            self.set_variable('ValveDataCodec', np.array(ValveDataCodec.version, dtype=np.int16))  # not used right now
        elif self.mode == 'r':
            versions = dict([(k, tuple(v)) for k, v in zip(['version', 'aquaduct_version'],
                                                           (self.get_variable('version'),
                                                            self.get_variable('aquaduct_version')))])
            check_versions(versions)

    def close(self):
        self.data_file.close()

    def get_variable(self, name, copy=True):
        # print name
        if copy:
            return self.data_file.variables[name][:].copy()
        return self.data_file.variables[name]

    def set_variable(self, name, value):
        # print name, value
        assert self.mode == "w"
        # value has to be ndarray
        assert isinstance(value, np.ndarray)
        # create dimensions
        dimensions = []
        for nr, d in enumerate(value.shape):
            dimensions.append("%s%d" % (name, nr))
            self.data_file.createDimension(dimensions[-1], d)
        # create variable
        if hasattr(netcdf, 'Dataset'):
            v = self.data_file.createVariable(name, value.dtype, tuple(dimensions), zlib=True)
        else:
            v = self.data_file.createVariable(name, value.dtype, tuple(dimensions))
        # fill variable
        v[:] = value

    def dump(self, **kwargs):
        for name, value in kwargs.items():
            for nname, vvalue in ValveDataCodec.encode(name, value):
                self.set_variable(nname, vvalue)

    @dictify
    def load(self):
        # names of data objects
        names = list(
            set([name.split('.')[0] for name in list(self.data_file.variables.keys()) if name not in self.not_variable]))
        all_names = [name for name in list(self.data_file.variables.keys()) if name not in self.not_variable]
        for name in names:
            # if 'inls' in name: continue
            # read all parts of object and decode
            this_object = dict(((n, self.get_variable(n, copy=False)) for n in all_names if
                                name == n or (name + '.') == n[:len(name) + 1]))
            yield name, ValveDataCodec.decode(name, this_object)
            for k in list(this_object.keys()):
                v = this_object.pop(k)[:]
                if hasattr(v, 'copy'):
                    v.copy()
                del v


################################################################################

class ValveDataAccess_numpy(ValveDataAccess):
    pass


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
                        for path_name, path in value.items():
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
        if False:  # hasattr(value,'simple_dump'):
            pickle.dump({name: value.simple_dump()}, self.data_file)
        else:
            pickle.dump({name: value}, self.data_file)
