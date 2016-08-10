import numpy as np

#from matplotlib.colors import colorConverter


from aqueduct.traj.paths import GenericPathTypeCodes as gptc
from aqueduct.traj.paths import PathTypesCodes as ptc
from aqueduct.visual.cmaps import default as default_cmap
from aqueduct.utils.helpers import zip_zip, is_number


_cl2rgba = {'b': (0.0, 0.0, 1.0, 1.0),
            'c': (0.0, 0.75, 0.75, 1.0),
            'g': (0.0, 0.5, 0.0, 1.0),
            'k': (0.0, 0.0, 0.0, 1.0),
            'm': (0.75, 0, 0.75, 1.0),
            'r': (1.0, 0.0, 0.0, 1.0),
            'y': (0.75, 0.75, 0, 1.0)}


#cc = lambda c, alpha=1.0: colorConverter.to_rgb(c)

def cc_safe(c):
    #color converter
    if c in 'rgbcmyk':
        c = _cl2rgba[c]
    c = tuple(c)
    assert len(c) in [3,4], 'Color should be given either as one letter (rgbcmyk) or as rgb or rgba vector.'
    for e in c:
        assert is_number(e), 'Color vector has to be specified as numbers.'
        assert e >=0 and e <= 1, 'Color vector elements have to be in range of 0 to 1.'
    return c[:3]

def cc(c):
    #color converter faster
    if c in 'rgbcmyk':
        c = _cl2rgba[c]
    return c[:3]



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

def get_cmap(size):
    return [e[0][0] for e in zip_zip(default_cmap,N=size)]

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
        self.cmap = get_cmap(self.cm_size)

    def __call__(self, node):
        if node > 0 and node <= self.size:
            return self.cmap[int(np.round(self.cm_size * f_like(node)))][:3]
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


