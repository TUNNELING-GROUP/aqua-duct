
# this modlue is a prototype and have to be revritten


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


part2type_dict = {0: GenericPathTypeCodes.scope_name,
                  1: GenericPathTypeCodes.object_name,
                  2: GenericPathTypeCodes.scope_name}

parts = (0, 1, 2)


def create_master_spath(spaths, smooth=None, resid=0, ctype=None, bias_long=5, heartbeat=None):

    def beat():
        if heartbeat is not None:
            heartbeat()

    # first check what is the size of paths in all parts and normalize and then scale them
    sizes = []
    for part in parts:
        # lengths of all paths of part part
        lens = np.array([float(len(sp.types[part])) for sp in spaths])
        if np.max(lens) > 0:
            lens /= np.max(lens) # normalization
            lens = lens ** bias_long # scale them by increasing weights of long paths
        if sum(lens) == 0:
            sizes.append(0)
        else:
            # weighted average by paths lengths
            sizes.append(int(np.average([len(sp.types[part]) for sp in spaths],0,lens)))
    full_size = sum(sizes) # total size (desired)

    beat() # touch progress bar

    # get total lenghts of all paths - needed as weights in averaging
    lens = np.array([float(sp.size) for sp in spaths])
    # if ctype in #:# and not 0 and not None then take object part only length
    if ctype is not None:
        if ctype.input is not None:
            if ctype.input > 0:
                if ctype.input == ctype.output:
                    lens = np.array([float(len(sp.types_object)) for sp in spaths])
    # normalize and incearse weight of long paths (depends ob ctype)
    if np.max(lens) > 0:
        lens /= np.max(lens) # normalize
        lens = lens ** bias_long # bias to long paths

    # containers for coords, types and widths of master path
    coords = []
    types = []
    widths = []
    # loop over zip zipped [smooth] coords of all paths and gtypes with size set to full_size
    for coords_zz, types_zz in zip(zip_zip(*[sp.get_coords_cont(smooth=smooth) for sp in spaths], N=full_size),
                                   zip_zip(*[sp.gtypes_cont for sp in spaths], N=full_size)):
        # make lens_zz which are lens corrected to the lenghts of coords_zz and normalized to zip_zip number of obejcts
        lens_zz = []
        for l, coord_z in zip(lens, coords_zz):
            if len(coord_z) > 0:
                lens_zz.append([float(l) / len(coord_z)] * len(coord_z)) # normalize and correct lengths
            else:
                # lens_zz.append([float(l)] * len(coord_z))
                lens_zz.append([])
        # concatenate zip_zip coords and lens
        coords_zz_cat = list(concatenate(*coords_zz))
        lens_zz_cat = list(concatenate(*lens_zz))
        # average coords_zz_cat using weights of lens_zz_cat
        coords.append(np.average(coords_zz_cat, 0, lens_zz_cat))
        # calculate widths
        if len(coords_zz) > 1:
            # try tu use weighted distance - wminkowski with p=2 is equivalent to weighted euclidean
            #id_of_max = np.argmax(pdist(coords_zz_cat, 'wminkowski', p=2, w=lens_zz_cat))
            #widths.append(pdist(coords_zz_cat, 'euclidean')[id_of_max])
            widths.append(np.mean(pdist(coords_zz_cat, 'euclidean')))
        else:
            widths.append(0.)
        # concatenate zip_zip gtypes
        types_zz_cat = list(concatenate(*types_zz))
        # pick correct type..., check distance of coords[-1] to coords_zz_cat
        types_cdist = cdist(np.matrix(coords[-1]),coords_zz_cat,metric='euclidean')
        types.append(types_zz_cat[np.argmin(types_cdist)])
        #types.append(decide_on_type(Counter(types_zz_cat)))

        beat()  # touch progress bar

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
    beat()  # touch progress bar
    sp = list(yield_single_paths([gp]))[0]
    beat()  # touch progress bar
    mp = MasterPath(sp)
    mp.add_width(widths)
    return mp




class MasterTrace(object):

    def __init__(self):
        pass
