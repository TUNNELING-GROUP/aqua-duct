import numpy as np
from scipy.spatial.distance import cdist, pdist ,squareform
from collections import Counter
from aqueduct.traj.paths import GenericPathTypeCodes, GenericPaths, yield_single_paths, MasterPath
from aqueduct.utils.helpers import list_blocks_to_slices
from aqueduct.geom import traces


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

def compress_zip(*args):
    ns = map(float,map(len,args))
    N = int(min(ns))
    position = [0.]*len(args)
    for n in range(N):
        this_yield = []
        next_position = [float(len(a))/N + p for a,p in zip(args,position)]
        for a,p,np in zip(args,position,next_position):
            if n + 1 == N:
                this_yield.append(a[int(p):])
            else:
                this_yield.append(a[int(p):int(np)])
        yield tuple(this_yield)
        position = next_position


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

def get_weights_(spaths,smooth=None):
    # max len
    max_len = [sp.get_length_cont(smooth=smooth)[-1] for sp in spaths]
    arg_max_len = np.argmax(max_len)
    max_len = max(max_len)
    # get weights as both lengths
    weights = np.array([sz for sz in strech_zip(*[sp.get_length_both_cont(smooth=smooth,normalize=max_len) for sp in spaths])])
    # add 0.5 - len_both of the lognest path
    length_both_of_max_len = 0.5 - spaths[arg_max_len].get_length_both_cont(smooth=smooth,normalize=max_len)

    weights  = np.array([lboml+w for w,lboml in strech_zip(weights,length_both_of_max_len)])
    return weights**10

def get_mean_coord_(coords,l):
    # l >> 0
    coord0 = np.median(coords,0)
    # l >> 1
    coord5 = coords[np.argmax(cdist(coords,np.matrix(coord0)))]

    return np.average([coord0,coord5],0,[1-l,l])

def concatenate(*args):
    for a in args:
        for e in a:
            yield e

def create_master_spath(spaths,smooth=None,resid=0,ctype=None):

    N = len(spaths)

    # coords:
    coords = []
    types = []
    for part in (0,1,2):
        #coords += [np.mean(list(concatenate(*cz)), 0) for cz in compress_zip(*[sp.get_coords(smooth=smooth)[part] for sp in spaths])]
        lens = np.array([float(len(sp.get_coords()[part])) for sp in spaths])
        if np.max(lens) > 0:
            lens /= np.max(lens)
            lens = lens**5
        coords_ = np.array([np.average(sz, 0,lens) for sz in strech_zip(*[sp.get_coords(smooth=smooth)[part] for sp in spaths])])
        #coords_ = traces.LinearizeRecursiveVector(0.05236)(coords_).tolist()
        coords.extend(coords_)
        types.extend(({0: 's', 1: 'c', 2: 's'}[part])*len(coords_)) # FIXME: magic constants!

    coords = np.array(coords)






    #weights = get_weights_(spaths,smooth=smooth)

    #coords = np.array([np.median(sz,0) for sz in strech_zip(*[sp.get_coords_cont(smooth=smooth) for sp in spaths])])
    #coords = np.array([np.average(sz,0,weights=w) for sz,w in zip(strech_zip(*[sp.get_coords_cont(smooth=smooth) for sp in spaths]),weights)])
    #max_len = [sp.get_length_cont(smooth=smooth)[-1] for sp in spaths]
    #lengths = np.array([np.max(sz) for sz in strech_zip(*[sp.get_length_both_cont(smooth=smooth,normalize=max_len) for sp in spaths])])



    #coords = np.array([get_mean_coord_(sz,l*2) for sz,l in zip(strech_zip(*[sp.get_coords_cont(smooth=smooth) for sp in spaths]),lengths)])


    # firs let's create types
    '''
    # types distribution
    td = [simple_types_distribution(sp.gtypes_cont) for sp in spaths]
    td = np.median(td,0)
    td.shape = (1,3)
    # get types as they are

    types = []
    for part in (0,1,2):
        #types += [Counter(concatenate(*cz)) for cz in compress_zip(*[sp.gtypes_cont[part] for sp in spaths])]
        types += [Counter(sz) for sz in strech_zip(*[sp.gtypes_cont[part] for sp in spaths])]

    def evaluate_types_counter(list_types_counter,s2o_treshold=0.5):
        resulted_types = map(lambda t: decide_on_type(t,s2o_treshold=s2o_treshold),list_types_counter)
        return np.array(simple_types_distribution(resulted_types))

    types_evaluation = np.array([evaluate_types_counter(types,s2o_treshold=s2o_t) for s2o_t in np.linspace(0,1,N+1)])

    s2o_treshold = np.linspace(0,1,N+1)[np.argmin(cdist(types_evaluation,td))]
    '''

    '''
    types = []
    for part in (0,1,2):
        types += [decide_on_type(Counter(concatenate(*cz)),s2o_treshold=s2o_treshold) for cz in compress_zip(*[sp.gtypes_cont for sp in spaths])]
    '''

    '''
    types = []
    for part in (0,1,2):
        types += [{0:'s',1:'c',2:'s'}[part] for sz in strech_zip(*[sp.gtypes[part] for sp in spaths])]
    '''


    #coords = np.array([np.mean(list(concatenate(*cz)), 0) for cz in compress_zip(*[sp.get_coords_cont(smooth=smooth) for sp in spaths])])



    if N == 1:
        width = np.zeros(len(coords))
    else:
        width = []
        for part in (0,1,2):
            #width += [np.median(pdist(list(concatenate(*cz)))) for cz in compress_zip(*[sp.get_coords(smooth=smooth)[part] for sp in spaths])]
            width += [np.median(pdist(sz)) for sz in strech_zip(*[sp.get_coords(smooth=smooth)[part] for sp in spaths])]
        width = np.array(width)


    #print coords

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
    #return gp
    sp = list(yield_single_paths([gp]))[0]
    mp = MasterPath(sp)
    mp.add_width(width)
    return mp




class MasterTrace(object):

    def __init__(self):
        pass

