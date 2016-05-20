import numpy as np
from scipy.spatial.distance import cdist, pdist ,squareform
from collections import Counter
from aqueduct.traj.paths import GenericPathTypeCodes, GenericPaths, yield_single_paths, MasterPath
from aqueduct.utils.helpers import list_blocks_to_slices

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




def decide_on_type(cont, s2o_treshold=0.5):
    # possible types are:
    #GenericPathTypeCodes.object_name
    #GenericPathTypeCodes.scope_name

    o = 0
    if cont.has_key(GenericPathTypeCodes.object_name):
        o = cont[GenericPathTypeCodes.object_name]
    s = 0
    if cont.has_key(GenericPathTypeCodes.scope_name):
        s = cont[GenericPathTypeCodes.scope_name]
    # get s to o value
    if s == o:
        s2o = 0.5
    else:
        if o == 0:
            s2o = 1.
        else:
            s2o = float(s)/(o+s)
    # decide on type
    #print s,o
    #print s2o, {True:'sco',False:'obj'}[s2o >= s2o_treshold],cont
    if s2o >= s2o_treshold and s > 0:
        return GenericPathTypeCodes.scope_name
    return GenericPathTypeCodes.object_name

def simple_types_distribution(types):
    # possible types are:
    #GenericPathTypeCodes.object_name
    #GenericPathTypeCodes.scope_name
    td_in, td_obj, td_out = 0,0,0
    sls = list(list_blocks_to_slices(types))
    if GenericPathTypeCodes.scope_name in types[sls[0]]:
        # this is input part
        td_in = len(types[sls[0]])
    if GenericPathTypeCodes.scope_name in types[sls[-1]]:
        # this is output part
        td_out = len(types[sls[-1]])
    # the rest is object
    td_obj = len(types) - td_in - td_out
    return map(lambda x: float(x)/len(types),(td_in,td_obj,td_out))

def create_master_spath(spaths,smooth=None,resid=0,ctype=None):

    N = len(spaths)

    coords = np.array([np.median(sz,0) for sz in strech_zip(*[sp.get_coords_cont(smooth=smooth) for sp in spaths])])
    if N == 1:
        width = np.zeros(len(coords))
    else:
        width = np.array([np.median(pdist(sz)) for sz in strech_zip(*[sp.get_coords_cont(smooth=smooth) for sp in spaths])])

    types_dist = [simple_types_distribution(sp.gtypes_cont) for sp in spaths]

    #print "types dist:"
    #print np.median(types_dist,0),sum(np.median(types_dist,0))
    #print np.mean(types_dist,0),sum(np.mean(types_dist,0))
    td = np.median(types_dist,0)
    td.shape = (1,3)

    # find optimal treshold
    types = [Counter(sz) for sz in strech_zip(*[sp.gtypes_cont for sp in spaths])]

    def evaluate_types_counter(list_types_counter,s2o_treshold=0.5):
        resulted_types = map(lambda t: decide_on_type(t,s2o_treshold=s2o_treshold),list_types_counter)
        return np.array(simple_types_distribution(resulted_types))

    types_evaluation = np.array([evaluate_types_counter(types,s2o_treshold=s2o_t) for s2o_t in np.linspace(0,1,N+1)])

    #print types_evaluation
    #print cdist(types_evaluation,td)
    #print np.argmin(cdist(types_evaluation,td))
    #print np.linspace(0,1,N+1)[np.argmin(cdist(types_evaluation,td))]
    s2o_treshold = np.linspace(0,1,N+1)[np.argmin(cdist(types_evaluation,td))]

    types = [decide_on_type(Counter(sz),s2o_treshold=s2o_treshold) for sz in strech_zip(*[sp.gtypes_cont for sp in spaths])]


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

