# -*- coding: utf-8 -*-
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import colorConverter

from aqueduct.geom import traces
from aqueduct.traj.paths import GenericPathTypeCodes as gptc
from aqueduct.traj.paths import PathTypesCodes as ptc
from aqueduct.utils.helpers import list_blocks_to_slices

# matplotlib specific
# mpl color converter
# cc = lambda c,alpha=1.0 : colorConverter.to_rgba(c,alpha=alpha)
cc = lambda c, alpha=1.0: colorConverter.to_rgb(c)


def get_cmap(name, size=None):
    return plt.get_cmap(name, lut=size)


class ColorMapDistMap(object):
    default_cm_size = 256

    grey = (0.5, 0.5, 0.5, 1)

    def __init__(self, name='hsv', size=None):
        # size is number of nodes to be maped to distinguistive colors
        self.size = size
        self.cm_size = self.default_cm_size
        while (self.cm_size < self.size):
            self.cm_size *= 1.1
            self.cm_size = int(np.ceil(self.cm_size))
        # get size
        self.cmap = get_cmap(name, self.cm_size)

    def __call__(self, node):
        if node > 0 and node <= self.size:
            return self.cmap(int(np.round(self.cm_size * f_like(node))))[:3]
        # return grey otherwise
        return self.grey[:3]


def f_like(n):
    if n == 1:
        return 0.0
    if n == 2:
        return 0.5
    n -= 1
    order = np.floor(np.log(n) / np.log(2))
    parts = 2 ** order
    current = n - parts
    return 0.5 / parts + 1. / parts * current


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


_dcc_is = ptc.path_in_code + gptc.scope_name
_dcc_cc = ptc.path_object_code + gptc.object_name
_dcc_cs = ptc.path_object_code + gptc.scope_name
_dcc_os = ptc.path_out_code + gptc.scope_name

_dcc_i = ptc.path_in_code
_dcc_c = ptc.path_object_code
_dcc_o = ptc.path_out_code

_default_color_codes = {_dcc_is: 'r',
                        _dcc_cc: 'g',
                        _dcc_cs: 'y',
                        _dcc_os: 'b',
                        _dcc_i: 'r',
                        _dcc_c: 'g',
                        _dcc_o: 'b'}

default_color_codes = _default_color_codes


def color_codes(code, custom_codes=None):
    if custom_codes is None:
        return default_color_codes[code]
    else:
        return custom_codes[code]


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
        raise NotImplementedError('This is base class.')

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
