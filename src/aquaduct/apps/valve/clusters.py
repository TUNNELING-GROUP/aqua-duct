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

from functools import wraps

import numpy as np
from scipy.stats import ttest_ind
from aquaduct.traj.inlets import InletTypeCodes
from aquaduct.apps.valve.spath import spath_lenght_total_info, spath_frames_info
from aquaduct.geom import hdr


def cluster_id_header():
    return ['Cluster'], ['%7d']


def add_cluster_id_head(gen):
    sph, splt = cluster_id_header()

    @wraps(gen)
    def patched(*args, **kwargs):
        h, lt = gen(*args, **kwargs)
        return sph + h, splt + lt

    return patched


def add_cluster_id(gen):
    @wraps(gen)
    def patched(cluster, something, add_id=True, *args, **kwargs):
        line = gen(cluster, something, *args, **kwargs)
        if add_id:
            line = [int(cluster)] + line
        return line

    return patched


@add_cluster_id_head
def clusters_inlets_header():
    header = 'Size SInp IInp IOut SOut'.split()
    header = 'Size INCOMING OUTGOING'.split()
    line_template = ['%7d'] + ['%8d'] * (len(header) - 1)
    return header, line_template


@add_cluster_id
def clusters_inlets(cluster, inlets):
    line = [inlets.size]
    line.append(inlets.lim2types([InletTypeCodes.surface_incoming]).size)
    # line.append(inlets.lim2types([InletTypeCodes.internal_incoming]).size)
    # line.append(inlets.lim2types([InletTypeCodes.internal_outgoing]).size)
    line.append(inlets.lim2types([InletTypeCodes.surface_outgoing]).size)
    return line


@add_cluster_id_head
def clusters_area_header():
    header = ('D' + ' D'.join(map(str, list(range(100, 85, -5)) + list(range(80, 40, -10))))).split()
    line_template = ['%8.2f'] * len(header)
    return header, line_template


@add_cluster_id
def clusters_area(cluster, inlets, points=10, expand_by=1):
    # HDR
    h = hdr.HDR(np.array(inlets.coords), points=points, expand_by=expand_by, center_of_system=inlets.center_of_system)
    line = list()
    for fraction in list(range(100, 85, -5)) + list(range(80, 40, -10)):
        line.append(h.area(fraction=fraction / 100.))
    return line


@add_cluster_id_head
def clusters_stats_prob_header():
    header = 'IN-OUT diff N IN-OUT_prob diff_prob N_prob'.split()
    line_template = ['%8d'] * int((len(header) / 2)) + ['%12.2f'] * int((len(header) / 2))
    # header += 'IN_len OUT_len Both_len'.split()
    # line_template += ['%9.1f'] * 3
    return header, line_template


@add_cluster_id
def clusters_stats_prob(cluster, sp_ct):
    # calculates probabilities of some events for cluster
    # X:X transition - io
    # X:? and ?:X transition - d
    # X:N and N:X transition - N
    line = []
    io, d, N = 0, 0, 0
    in_len, out_len, tot_len = 0., 0., 0.
    in_n, out_n, tot_n = 0, 0, 0
    for sp, ct in sp_ct:
        ct = ct.generic.clusters
        assert cluster in ct
        if cluster == ct[0] and cluster == ct[1]:
            io += 1
        elif None in ct:
            N += 1
        else:
            d += 1
    line += [io, d, N]
    summa = float(sum([io, d, N]))
    line += [x / summa if summa else float('nan') for x in [io, d, N]]
    return line


@add_cluster_id_head
def clusters_stats_len_header():
    header = 'X->Obj Obj->X p-value X->ObjMin X->ObjMinID Obj->XMin Obj->XMinID'.split()
    line_template = (['%9.1f'] * 2) + ['%9.4f'] + ['%9.1f', '%11s'] * 2
    return header, line_template


@add_cluster_id
def clusters_stats_len(cluster, sp_ct):
    line = []
    in_len, out_len = [], []
    in_len_min, out_len_min = float('inf'), float('inf')
    in_len_min_id, out_len_min_id = None, None
    for sp, ct in sp_ct:
        ct = ct.generic.clusters
        assert cluster in ct
        lens = spath_lenght_total_info(sp, add_id=False, total=False, totalonly=False)
        # tot,in,obj,out
        if cluster == ct[0]:
            if not np.isnan(lens[0]):
                in_len.append(lens[0])
                if lens[0] < in_len_min:
                    in_len_min = lens[0]
                    in_len_min_id = str(sp.id)
        if cluster == ct[1]:
            if not np.isnan(lens[-1]):
                out_len.append(lens[-1])
                if lens[-1] < out_len_min:
                    out_len_min = lens[-1]
                    out_len_min_id = str(sp.id)

    if len(in_len):
        line.append(np.mean(in_len))
    else:
        line.append(float('nan'))
    if len(out_len):
        line.append(np.mean(out_len))
    else:
        line.append(float('nan'))
    if len(in_len) > 1 and len(out_len) > 1:
        line.append(ttest_ind(in_len, out_len)[-1])  # this is supposed to return p-value
    else:
        line.append(float('nan'))

    line += [in_len_min, in_len_min_id, out_len_min, out_len_min_id]

    return line


@add_cluster_id_head
def clusters_stats_steps_header():
    header = 'X->Obj Obj->X p-value X->ObjMin X->ObjMinID Obj->XMin Obj->XMinID'.split()
    line_template = (['%9.1f'] * 2) + ['%9.4f'] + ['%9.1f', '%11s'] * 2
    return header, line_template


@add_cluster_id
def clusters_stats_steps(cluster, sp_ct):
    line = []
    in_len, out_len = [], []
    in_len_min, out_len_min = float('inf'), float('inf')
    in_len_min_id, out_len_min_id = None, None
    for sp, ct in sp_ct:
        ct = ct.generic.clusters
        assert cluster in ct
        lens = spath_frames_info(sp, add_id=False, total=False)
        # tot,in,obj,out
        if cluster == ct[0]:
            if not np.isnan(lens[0]):
                in_len.append(lens[0])
                if lens[0] < in_len_min:
                    in_len_min = lens[0]
                    in_len_min_id = str(sp.id)
        if cluster == ct[1]:
            if not np.isnan(lens[-1]):
                out_len.append(lens[-1])
                if lens[-1] < out_len_min:
                    out_len_min = lens[-1]
                    out_len_min_id = str(sp.id)

    if len(in_len):
        line.append(np.mean(in_len))
    else:
        line.append(float('nan'))
    if len(out_len):
        line.append(np.mean(out_len))
    else:
        line.append(float('nan'))
    if len(in_len) > 1 and len(out_len) > 1:
        line.append(ttest_ind(in_len, out_len)[-1])  # this is supposed to return p-value
    else:
        line.append(float('nan'))

    line += [in_len_min, in_len_min_id, out_len_min, out_len_min_id]

    return line
