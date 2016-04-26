import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import colorConverter
from mpl_toolkits.mplot3d import Axes3D

from aqueduct.geom import traces
from aqueduct.utils.helpers import list_blocks_to_slices

# matplotlib specific
# mpl color converter
# cc = lambda c,alpha=1.0 : colorConverter.to_rgba(c,alpha=alpha)
cc = lambda c, alpha=1.0: colorConverter.to_rgb(c)


def get_cmap(name, size=None):
    return plt.get_cmap(name, lut=size)


class ColorMapDistMap(object):
    dist = 1

    grey = (0.5, 0.5, 0.5, 1)

    def __init__(self, name='hsv', size=None):
        # size is number of nodes to be maped to distinguistive colors
        self.size = size
        # lut should be appropriately bigger - 10 times
        lut = size * self.dist + 1
        # get size
        self.cmap = get_cmap(name, lut)

    def __call__(self, node):
        if node > 0 and node <= self.size:
            # get color
            return self.cmap(int(node) * self.dist)[:3]
        # return grey otherwise
        return self.grey[:3]


def yield_spath_len_and_smooth_diff_in_types_slices(sp, smooth=None, smooth_len=None, smooth_diff=None):
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
    #ldif = np.array(range(len(dif))) # same OX

    etypes = sp.etypes_cont

    for sl in list_blocks_to_slices(etypes):
        etype = etypes[sl]
        ld = ldif[sl]
        sd = dif[sl]
        while len(etype) > len(ld):
            etype.pop(-1)
        yield ld, sd, etype




default_color_codes = {'is': 'r',
                       'cc': 'g',
                       'cs': 'y',
                       'os': 'b'}


def color_codes(code, custom_codes=None):
    if custom_codes is None:
        return default_color_codes[code]
    else:
        return custom_codes[code]


def plot_spath_spectrum(sp, **kwargs):
    lsdt = list(yield_spath_len_and_smooth_diff_in_types_slices(sp, **kwargs))

    n = len(lsdt)

    last_l = None
    last_sd = None
    last_t = None
    last_color = None
    for nr, (l, sd, t) in enumerate(lsdt):
        color = color_codes(t[-1])
        if nr == 0:
            plt.plot(l, sd, color=color)
            last_l = l[-1]
            last_sd = sd[-1]
            last_color = color
        else:
            mid_l = (last_l + l[0]) / 2
            mid_sd = (last_sd + sd[0]) / 2
            if nr == n - 1:
                plt.plot([last_l, mid_l], [last_sd, mid_sd], color=last_color)
                plt.plot([mid_l, l[0]], [mid_sd, sd[0]], color=color)
                plt.plot(l, sd, color=color)
            else:
                plt.plot([last_l, mid_l], [last_sd, mid_sd], color=last_color)
                plt.plot([mid_l, l[0]], [mid_sd, sd[0]], color=color)
                plt.plot(l, sd, color=color)
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
