# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2016-2017  Tomasz Magdziarz, Alicja Płuciennik, Michał Stolarczyk <info@aquaduct.pl>
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

import numpy as np

import matplotlib.pyplot as plt

from aquaduct.geom import traces
from aquaduct.utils.helpers import list_blocks_to_slices
from aquaduct.visual.helpers import color_codes, cc


# matplotlib specific
# mpl color converter
# cc = lambda c,alpha=1.0 : colorConverter.to_rgba(c,alpha=alpha)


def yield_spath_len_and_smooth_diff_in_types_slices(sp, smooth=None, smooth_len=None, smooth_diff=None, types='etypes'):
    if smooth is not None:
        smooth_len = smooth
        smooth_diff = smooth
    # len
    coords = sp.get_coords_cont(smooth=smooth_len)
    dif = traces.diff(coords)
    ldif = np.cumsum(dif)
    # diff
    if smooth_diff != smooth_len:
        coords = sp.get_coords_cont(smooth=smooth_diff)
        dif = traces.diff(coords)
    # ldif = np.array(range(len(dif))) # same OX

    if types == 'etypes':
        etypes = sp.etypes_cont
    elif types == 'types':
        etypes = sp.types_cont

    for sl in list_blocks_to_slices(etypes):
        etype = etypes[sl]
        ld = ldif[sl]
        sd = dif[sl]
        while len(etype) > len(ld):
            etype.pop(-1)
        yield ld, sd, etype


def plot_colorful_lines(x, y, c, **kwargs):
    sls = list(list_blocks_to_slices(c))
    n = len(sls)

    for nr, sl in enumerate(sls):
        a = x[sl]
        b = y[sl]
        color = c[sl][-1]
        if nr == 0:
            plt.plot(a, b, color=color, **kwargs)
            last_a = a[-1]
            last_b = b[-1]
            last_color = color
        else:
            mid_a = (last_a + a[0]) / 2
            mid_b = (last_b + b[0]) / 2
            if nr == n - 1:
                plt.plot([last_a, mid_a], [last_b, mid_b], color=last_color, **kwargs)
                plt.plot([mid_a, a[0]], [mid_b, b[0]], color=color, **kwargs)
                plt.plot(a, b, color=color, **kwargs)
            else:
                plt.plot([last_a, mid_a], [last_b, mid_b], color=last_color, **kwargs)
                plt.plot([mid_a, a[0]], [mid_b, b[0]], color=color, **kwargs)
                plt.plot(a, b, color=color, **kwargs)
                last_a = a[-1]
                last_b = b[-1]
                last_color = color


def spaths_spectra(spaths, **kwargs):
    spectra = []

    minx = None
    maxx = None
    for sp in spaths:
        xyc = list(spath_spectrum(sp, **kwargs))
        xy = np.array([(x, y) for x, y, c in xyc])
        c = list([c for x, y, c in xyc])
        x = xy[:, 0]
        y = xy[:, 1]
        spectra.append((y, c))
        if maxx is None:
            maxx = max(x)
        else:
            maxx = max(maxx, max(x))
        if minx is None:
            minx = min(x)
        else:
            minx = min(minx, min(x))
    for (y, c) in spectra:
        plot_colorful_lines(np.linspace(minx, maxx, len(c)), y, c)


def plot_spath_spectrum(sp, **kwargs):
    xyc = list(spath_spectrum(sp, **kwargs))
    xy = np.array([(x, y) for x, y, c in xyc])
    c = list([c for x, y, c in xyc])
    plot_colorful_lines(xy[:, 0], xy[:, 1], c)


def spath_spectrum(sp, **kwargs):
    lsdt = list(yield_spath_len_and_smooth_diff_in_types_slices(sp, **kwargs))

    n = len(lsdt)

    last_l = None
    last_sd = None
    last_t = None
    last_color = None
    for nr, (l, sd, t) in enumerate(lsdt):
        color = color_codes(t[-1])
        if nr == 0:
            for ll, ssdd in zip(l, sd):
                yield ll, ssdd, color
            # plt.plot(l, sd, color=color)
            last_l = l[-1]
            last_sd = sd[-1]
            last_color = color
        else:
            mid_l = (last_l + l[0]) / 2
            mid_sd = (last_sd + sd[0]) / 2
            if nr == n - 1:
                # for ll, ssdd in zip([last_l, mid_l], [last_sd, mid_sd]):
                #    yield ll, ssdd, last_color
                # plt.plot([last_l, mid_l], [last_sd, mid_sd], color=last_color)
                # for ll, ssdd in zip([mid_l, l[0]], [mid_sd, sd[0]]):
                #    yield ll, ssdd, color
                # plt.plot([mid_l, l[0]], [mid_sd, sd[0]], color=color)
                for ll, ssdd in zip(l, sd):
                    yield ll, ssdd, color
                    # plt.plot(l, sd, color=color)
            else:
                # for ll, ssdd in zip([last_l, mid_l], [last_sd, mid_sd]):
                #    yield ll, ssdd, last_color
                # plt.plot([last_l, mid_l], [last_sd, mid_sd], color=last_color)
                # for ll, ssdd in zip([mid_l, l[0]], [mid_sd, sd[0]]):
                #    yield ll, ssdd, color
                # plt.plot([mid_l, l[0]], [mid_sd, sd[0]], color=color)
                for ll, ssdd in zip(l, sd):
                    yield ll, ssdd, color
                # plt.plot(l, sd, color=color)
                last_l = l[-1]
                last_sd = sd[-1]
                last_color = color


def showit(gen):
    def patched(*args, **kwargs):
        obj = gen(*args, **kwargs)
        plt.show()
        return obj

    return patched


def get_ax3d(fig, sub=111):
    return fig.add_subplot(sub, projection='3d')


class SimpleTracePlotter(object):
    def plot_line(self, coords, color, **kwargs):
        raise NotImplementedError('This is abstract class. Missing implementaion in a child class.')

    '''
    def path_trace(self,path,color=('r','g','b'),
                   plot_in=True,
                   plot_object=True,
                   plot_out=True,
                   **kwargs):
        raise NotImplementedError('This is base class.')
    '''

    def single_trace(self, coords, color='r', **kwargs):
        # coords is a trace
        # color is a single color
        color = cc(color)
        coords = np.array(coords)
        # call plot_line
        self.plot_line(coords, color, **kwargs)

    def path_trace(self, path, color=('r', 'g', 'b'),
                   plot_in=True,
                   plot_object=True,
                   plot_out=True,
                   **kwargs):
        # path is a tuple of length 3, its elements represent in,object, out parts of path
        # color is a tuple of length 3, its elements correspond to colors of consecutive path parts
        color = map(cc, color)
        for nr, trace in enumerate(traces.midpoints(path)):
            # mid points!
            if len(trace) > 0:
                if (nr == 0 and plot_in) or (nr == 1 and plot_object) or (nr == 2 and plot_out):
                    self.single_trace(trace, color=color[nr], **kwargs)


class SimpleProteinPlotter(SimpleTracePlotter):
    def protein_trace(self, protein, smooth=None, color=('c', 'm', 'y'), **kwargs):
        # assumes protein is reader object
        # TODO: iterate over chains?
        bb = protein.parse_selection("protein and backbone")
        coords = bb.get_positions()
        cdiff = traces.diff(coords)
        split = np.argwhere(cdiff > 2.5)  # TODO: magic constant
        ns = len(split)
        if ns == 0:
            self.single_trace(smooth(coords), color=color[0], **kwargs)  # TODO: color conversion is buggy
        else:
            split.shape = (ns,)
            for nr, csplit in enumerate([0] + split.tolist()):
                cc = color[nr % len(color)]
                if nr == ns:
                    scoords = coords[csplit + 1:]
                elif nr == 0:
                    scoords = coords[csplit:split[nr]]
                else:
                    scoords = coords[csplit + 1:split[nr]]
                if smooth:
                    scoords = smooth(scoords)
                self.single_trace(scoords, color=color, **kwargs)


class SimplePathPlotter(SimpleTracePlotter):
    def single_path_traces(self, spaths, smooth=None, color=('r', 'g', 'b'), **kwargs):
        for spath in spaths:
            self.path_trace(spath.get_coords(smooth), color=color, **kwargs)


class MPLTracePlotter(SimplePathPlotter, SimpleProteinPlotter):
    @showit
    def init_ax(self):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection='3d')

        self.fig.subplots_adjust(left=0, bottom=0, right=1, top=1)
        self.fig.set_facecolor('w')

        self.ax.set_axis_bgcolor('none')
        self.ax.axis('off')

    @showit
    def plot_line(self, coords, color, **kwargs):
        self.ax.plot3D(coords[:, 0],
                       coords[:, 1],
                       coords[:, 2],
                       c=color, **kwargs)

    @showit
    def scatter(self, coords, **kwargs):
        self.ax.scatter3D(coords[:, 0],
                          coords[:, 1],
                          coords[:, 2],
                          **kwargs)
