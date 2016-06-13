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

def zip_zip(*args,**kwargs):
    ns = map(float,map(len,args))
    if 'N' in kwargs.keys():
        N = kwargs['N']
    else:
        N = int(min(ns))
    position = [0.]*len(args)
    for n in range(N):
        this_yield = []
        next_position = [float(len(a))/N + p for a,p in zip(args,position)]
        for a,p,np in zip(args,position,next_position):
            ip = int(p)
            inp = int(np)
            if n + 1 == N:
                this_yield.append(a[ip:])
            else:
                if ip == inp:
                    inp += 1
                this_yield.append(a[ip:inp])
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
    max_len = [sp.get_distance_cont(smooth=smooth)[-1] for sp in spaths]
    arg_max_len = np.argmax(max_len)
    max_len = max(max_len)
    # get weights as both lengths
    weights = np.array([sz for sz in strech_zip(*[sp.get_distance_both_cont(smooth=smooth, normalize=max_len) for sp in spaths])])
    # add 0.5 - len_both of the lognest path
    length_both_of_max_len = 0.5 - spaths[arg_max_len].get_distance_both_cont(smooth=smooth, normalize=max_len)

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

def create_master_spath(spaths, smooth=None, resid=0, ctype=None, bias_long=5, heartbeat=None):

    part2type_dict = {0: GenericPathTypeCodes.scope_name,
                      1: GenericPathTypeCodes.object_name,
                      2: GenericPathTypeCodes.scope_name}
    parts = (0, 1, 2)
    # first check what is the size of paths in all parts
    sizes = []
    for part in parts:
        lens = np.array([float(len(sp.types[part])) for sp in spaths])
        if np.max(lens) > 0:
            lens /= np.max(lens)
            lens = lens ** bias_long
        if sum(lens) == 0:
            sizes.append(0)
        else:
            sizes.append(int(np.average([len(sp.types[part]) for sp in spaths],0,lens)))
    #print ctype, sizes
    # now let's create coords, types and widths
    coords = []
    types = []
    widths = []
    for part in parts:
        #coords.append([np.mean(list(concatenate(*cz)), 0) for cz in zip_zip(*[sp.get_coords(smooth=smooth)[part] for sp in spaths],N=sizes[part])])

        lens = np.array([float(len(sp.types[part])) for sp in spaths])
        if np.max(lens) > 0:
            lens /= np.max(lens)
            lens = lens ** bias_long
        coords_ = []
        widths_ = []
        if sum(lens) != 0:
            for coords_zz in zip_zip(*[sp.get_coords(smooth=smooth)[part] for sp in spaths], N=sizes[part]):
                # make lens_zz which are lens corrected to the lenght of coord_z
                lens_zz = []
                for l,coord_z in zip(lens,coords_zz):
                    if len(coord_z) > 0:
                        lens_zz.append([float(l)/len(coord_z)] * len(coord_z))
                    else:
                        lens_zz.append([float(l)] * len(coord_z))
                #lens_zz = [[float(l)/len(coord_z)] * len(coord_z) for l,coord_z in zip(lens,coords_zz)]
                coords_zz_cat = list(concatenate(*coords_zz))
                lens_zz = list(concatenate(*lens_zz))
                coords_.append(np.average(coords_zz_cat, 0, lens_zz))
                if len(coords_zz) > 1:
                    widths_.append(np.max(pdist(coords_zz_cat, 'minkowski', p=2, w=lens_zz)))
                else:
                    widths_.append(0.)
        if heartbeat is not None:
            heartbeat()

        coords.append(coords_)
        widths.append(widths_)

        types.append([(part2type_dict[part])] * len(coords[-1]))

    # concatenate
    coords = list(concatenate(*coords))
    types = list(concatenate(*types))
    widths = list(concatenate(*widths))

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
    mp.add_width(widths)
    return mp




class MasterTrace(object):

    def __init__(self):
        pass

