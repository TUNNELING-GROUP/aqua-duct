import numpy as np
from scipy.spatial.distance import cdist, pdist ,squareform
from collections import Counter
from aqueduct.traj.paths import GenericPathTypeCodes, GenericPaths, yield_single_paths, MasterPath


def fit_trace_to_points(trace,points):
    dist = cdist(trace,points)
    points_min = np.argmin(dist,1).tolist()
    points_used = []
    for pm in points_min:
        if pm in points_used:
            continue
        points_used.append(pm)
    return points_used

def strech_zip(*args):

    ns = map(float,map(len,args))
    N = int(max(ns))
    for n in range(N):
        yield tuple([args[nr][int(cN/N*n)] for nr, cN in enumerate(ns)])




def decide_on_type(cont):
    # possible types are:
    #GenericPathTypeCodes.object_name
    #GenericPathTypeCodes.scope_name

    o = 0
    if cont.has_key(GenericPathTypeCodes.object_name):
        o = cont[GenericPathTypeCodes.object_name]
    s = 0
    if cont.has_key(GenericPathTypeCodes.scope_name):
        s = cont[GenericPathTypeCodes.scope_name]
    if s == o:
        return GenericPathTypeCodes.object_name
    return cont.most_common(1)[0][0]


def create_master_spath(spaths,smooth=None,resid=0,ctype=None):

    coords = np.array([np.median(sz,0) for sz in strech_zip(*[sp.get_coords_cont(smooth=smooth) for sp in spaths])])
    width = np.array([np.median(pdist(sz)) for sz in strech_zip(*[sp.get_coords_cont(smooth=smooth) for sp in spaths])])
    types = [decide_on_type(Counter(sz)) for sz in strech_zip(*[sp.gtypes_cont for sp in spaths])]
    frames = range(len(coords))

    min_pf = 0
    max_pf = len(coords) - 1

    if ctype is None:
        min_pf = None
        max_pf = None
    else:
        if ctype.input is not None:
            min_pf = None
        if ctype.output is not None:
            max_pf = None

    gp = GenericPaths(resid, min_pf=min_pf, max_pf=max_pf)
    for c,t,f in zip(coords,types,frames): # TODO: remove loop
        gp.add_type(f,t)
        gp.add_coord(c)
    sp = list(yield_single_paths([gp]))[0]
    mp = MasterPath(sp)
    mp.add_width(width)
    return mp




class MasterTrace(object):

    def __init__(self):
        pass

