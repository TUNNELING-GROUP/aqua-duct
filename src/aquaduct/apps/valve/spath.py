# -*- coding: utf-8 -*-

# This program is distributed in the hope that it will be useful,
# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2016-2018  Tomasz Magdziarz, Alicja Płuciennik, Michał Stolarczyk <info@aquaduct.pl>
# Copyright (C) 2019  Tomasz Magdziarz, Michał Banas <info@aquaduct.pl>
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


# Module contains mostly stuff need for stage V

from functools import wraps

import numpy as np
from aquaduct.geom import traces
from aquaduct.traj.paths import PassingPath


class SpathIdHeader(object):
    name = 'ID'
    format = '%9s'

    def __call__(self):
        return [self.name], [self.format]


'''
def spath_id_header():
    return ['ID'], ['%9s']
'''
spath_id_header = SpathIdHeader()


def spath_name_header():
    return ['RES'], ['%4s']


def add_path_id_head(gen):
    @wraps(gen)
    def patched(*args, **kwargs):
        sph, splt = list(zip(spath_id_header(), spath_name_header()))
        sph = [e[0] for e in sph]
        splt = [e[0] for e in splt]

        add_id = True
        if 'add_id' in kwargs:
            add_id = kwargs.pop('add_id')
        h, lt = gen(*args, **kwargs)
        if add_id:
            return sph + h, splt + lt
        return h, lt

    return patched


def add_path_id(gen):
    @wraps(gen)
    def patched(spath, *args, **kwargs):
        add_id = True
        if 'add_id' in kwargs:
            add_id = kwargs.pop('add_id')
        line = gen(spath, *args, **kwargs)
        if add_id:
            line = [spath.id, spath.id.name] + line
        return line

    return patched


def size_header():
    return ['Size', 'Size%'], ['%7d', '%6.2f']


def add_size_head(gen):
    sph, splt = size_header()

    @wraps(gen)
    def patched(*args, **kwargs):
        h, lt = gen(*args, **kwargs)
        return sph + h, splt + lt

    return patched


def add_size(gen):
    @wraps(gen)
    def patched(spaths, add_size=True, add_size_p100=None, *args, **kwargs):
        line = gen(spaths, *args, **kwargs)
        if add_size_p100 is not None:
            line = [len(spaths) / float(add_size_p100) * 100] + line
        if add_size:
            line = [len(spaths)] + line
        return line

    return patched


def ctype_id_header():
    return ['CType'], ['%7s']


def add_ctype_id_head(gen):
    sph, splt = ctype_id_header()

    @wraps(gen)
    def patched(*args, **kwargs):
        h, lt = gen(*args, **kwargs)
        return sph + h, splt + lt

    return patched


def add_ctype_id(gen):
    @wraps(gen)
    def patched(ctype, something, add_id=True, *args, **kwargs):
        line = gen(ctype, something, *args, **kwargs)
        if add_id:
            line = [str(ctype)] + line
        return line

    return patched


@add_path_id_head
def spath_basic_info_header():
    header = 'BeginF InpF ObjF ObjFS OutF EndF'.split()
    line_template = ['%7d'] * len(header)
    return header, line_template


@add_path_id
def spath_basic_info(spath):
    line = [spath.begins]
    if not isinstance(spath, PassingPath):
        line.extend(list(map(len, (spath.path_in, spath.path_object))))
        line.append(spath.path_object_strict_len)
        line.append(len(spath.path_out))
        # line.extend(map(len, (spath.path_in, spath.path_object, spath.path_out)))
    else:
        line += [float('nan')] * 4
    line.append(spath.ends)
    return line


@add_path_id_head
def spath_lenght_total_info_header(total=None):
    header = 'InpL ObjL OutL'.split()
    if total:
        header = ['TotL'] + header
    line_template = ['%9.1f'] * len(header)
    return header, line_template


@add_path_id
def spath_lenght_total_info(spath, totalonly=False, total=False):
    # total and totalonly are internal flags
    # to calculate: total len, in len, obj len, out len call:
    # total=True (this will calculate total)
    # to skip calculation of total call:
    # total=False, totalonly=False (this will calculate in, obj, and out lens)
    # to calculate total len call:
    # total=False, totalonly=True
    line = []
    if not total:
        if not totalonly:
            for t in traces.midpoints(spath.coords):
                if len(t) > 1:  # traces.length_step_std requires at least 2 points
                    line.append(traces.length_step_std(t)[0])
                else:
                    line.append(float('nan'))
            return line
        else:
            for t in traces.midpoints((spath.coords_cont,)):
                if len(t) > 1:  # traces.length_step_std requires at least 2 points
                    line.append(traces.length_step_std(t)[0])
                else:
                    line.append(float('nan'))
            return line
    line += spath_lenght_total_info(spath, add_id=False, total=False, totalonly=True)
    if not isinstance(spath, PassingPath):
        line += spath_lenght_total_info(spath, add_id=False, total=False, totalonly=False)
    else:
        line += [float('nan')] * 3
    return line


@add_path_id
def spath_frames_total_info(spath, totalonly=False, total=False):
    # total and totalonly are internal flags
    # to calculate: total len, in len, obj len, out len call:
    # total=True (this will calculate total)
    # to skip calculation of total call:
    # total=False, totalonly=False (this will calculate in, obj, and out lens)
    # to calculate total len call:
    # total=False, totalonly=True
    line = []
    if not total:
        if not totalonly:
            for t in spath.coords:
                line.append(len(t))
            return line
        else:
            for t in (spath.coords_cont,):
                line.append(len(t))
            return line
    line += spath_frames_total_info(spath, add_id=False, total=False, totalonly=True)
    if not isinstance(spath, PassingPath):
        line += spath_frames_total_info(spath, add_id=False, total=False, totalonly=False)
    else:
        line += [float('nan')] * 3
    return line


@add_path_id_head
def spath_steps_info_header(total=None):
    header = 'InpS InpStdS ObjS ObjStdS OutS OutStdS'.split()
    if total:
        header = ['TotS', 'TotStdS'] + header
    line_template = ['%8.2f', '%8.3f'] * int(len(header) / 2)
    return header, line_template


@add_path_id
def spath_steps_info(spath, total=None):
    line = []
    if not total:
        if not isinstance(spath, PassingPath):
            for t in traces.midpoints(spath.coords):
                if len(t) > 0:
                    line.extend(traces.length_step_std(t)[1:])
                else:
                    line.extend([float('nan'), float('nan')])
        else:
            line += [float('nan'), float('nan')] * 3
        return line
    line += spath_steps_info(spath, add_id=False, total=False)
    t = next(traces.midpoints((spath.coords_cont,)))
    if len(t) > 0:
        line = list(traces.length_step_std(t)[1:]) + line
    else:
        line = [float('nan'), float('nan')] + line
    return line


@add_path_id
def spath_frames_info(spath, total=None):
    line = []
    if not total:
        for t in spath.coords:
            line.append(len(t))
        return line
    line += spath_steps_info(spath, add_id=False, total=False)
    if not isinstance(spath, PassingPath):
        t = spath.coords_cont
        line = [len(t)] + line
    else:
        line = [float('nan')] + line
    return line


@add_path_id_head
def spath_ctype_header():
    header, line_template = ctype_id_header()
    return header, line_template


@add_path_id
def spath_ctype(spath, ctype=None):
    line = [str(ctype)]
    return line


@add_path_id_head
def spath_full_info_header(total=None):
    header = []
    line_template = []
    for h, lt in (spath_basic_info_header(add_id=False),
                  spath_lenght_total_info_header(add_id=False, total=total),
                  spath_steps_info_header(add_id=False, total=total),
                  spath_ctype_header(add_id=False)):
        header += h
        line_template += lt
    return header, line_template


@add_path_id
def spath_full_info(spath, ctype=None, total=None):
    line = []
    for l in (spath_basic_info(spath, add_id=False),
              spath_lenght_total_info(spath, add_id=False, total=total),
              spath_steps_info(spath, add_id=False, total=total),
              spath_ctype(spath, ctype=ctype, add_id=False)):
        line += l
    return line


@add_size_head
def spaths_lenght_total_header():
    header = 'Tot TotStd Inp InpStd Obj ObjStd Out OutStd'.split()
    line_template = ['%9.1f', '%9.2f'] * int(len(header) / 2)
    return header, line_template


@add_size
def spaths_length_total(spaths):
    line = []
    d4s = []
    for sp in spaths:
        if not isinstance(sp, PassingPath):
            d4s.append(spath_lenght_total_info(sp, add_id=False))
    d4s = np.array(d4s)
    if d4s.size:
        line.extend(np.mean(d4s, 0))
        line.extend(np.std(d4s, 0))
    else:
        line = [float('nan')] * 6
    # total
    d4s = []
    for sp in spaths:
        d4s.append(spath_lenght_total_info(sp, totalonly=True, add_id=False))
    d4s = np.array(d4s)
    if d4s.size:
        line.extend(np.mean(d4s, 0))
        line.extend(np.std(d4s, 0))
    else:
        line += [float('nan')] * 4

    return [line[6], line[7], line[0], line[3], line[1], line[4], line[2], line[5]]


@add_size
def spaths_frames_total(spaths):
    line = []
    d4s = []
    for sp in spaths:
        if not isinstance(sp, PassingPath):
            d4s.append(spath_frames_total_info(sp, add_id=False))
    d4s = np.array(d4s)
    if d4s.size:
        line.extend(np.mean(d4s, 0))
        line.extend(np.std(d4s, 0))
    else:
        line = [float('nan')] * 6
    # total
    d4s = []
    for sp in spaths:
        d4s.append(spath_frames_total_info(sp, totalonly=True, add_id=False))
    d4s = np.array(d4s)
    if d4s.size:
        line.extend(np.mean(d4s, 0))
        line.extend(np.std(d4s, 0))
    else:
        line += [float('nan')] * 4

    return [line[6], line[7], line[0], line[3], line[1], line[4], line[2], line[5]]


@add_ctype_id_head
def ctypes_spaths_info_header():
    header, line_template = spaths_lenght_total_header()
    return header, line_template


@add_ctype_id
def ctypes_spaths_info(ctype, spaths, show='len', add_size_p100=None):
    # show could be len or frames
    line = []
    if show == 'len':
        line += spaths_length_total(spaths, add_size_p100=add_size_p100)
    if show == 'frames':
        line += spaths_frames_total(spaths, add_size_p100=add_size_p100)
    return line


def plot_spaths_traces(spaths, spp=None, name=None, split=False, states=False, separate=False, smooth=None):
    if states or separate:
        spaths_iter = spaths
    else:
        spaths_iter = [spaths]
    state = None
    name_separate = ''
    for nr, sp in enumerate(spaths_iter):
        if states:
            state = nr + 1
        if separate:
            name_separate = '_%d' % (nr + 1)
        if states or separate:
            sp = [sp]
        if split:
            # FIXME: if path(s) is empty this will probably produce CGO_END without CGO_BEGIN
            spp.paths_trace(sp, name=name + '_in' + name_separate, plot_walk=False, plot_object=False, plot_out=False,
                            state=state, smooth=smooth)
            spp.paths_trace(sp, name=name + '_obj' + name_separate, plot_walk=False, plot_in=False, plot_out=False,
                            state=state, smooth=smooth)
            spp.paths_trace(sp, name=name + '_out' + name_separate, plot_walk=False, plot_in=False, plot_object=False,
                            state=state, smooth=smooth)
            spp.paths_trace(sp, name=name + '_walk' + name_separate, plot_in=False, plot_object=False, plot_out=False,
                            state=state, smooth=smooth)
        else:
            spp.paths_trace(sp, name=name + name_separate, state=state, smooth=smooth)


def plot_spaths_inlets(spaths, spp=None, name=None, states=False, separate=False, smooth=None):
    if states or separate:
        spaths_iter = spaths
    else:
        spaths_iter = [spaths]
    state = None
    name_separate = ''
    for nr, sp in enumerate(spaths_iter):
        if states:
            state = nr + 1
        if separate:
            name_separate = '_%d' % (nr + 1)
        if states or separate:
            sp = [sp]
        spp.paths_inlets(sp, name=name + name_separate, state=state, smooth=smooth)
