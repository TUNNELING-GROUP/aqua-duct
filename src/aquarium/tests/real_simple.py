
from aquarium import tests as aqtests
from aquarium.traj.reader import ReadAmberNetCDFviaMDA
from aquarium.traj.paths import GenericPath
from aquarium.utils import log

import multiprocessing as mp
import copy

def CHullCheck(point):
    return CHullCheck.chull.point_within(point)

def CHullCheck_init(args):
    CHullCheck.chull = copy.deepcopy(args[0])

def CHullCheck_pool(chull,threads=3):
    return mp.Pool(threads, CHullCheck_init, [(chull,)])

def CHullCheck_exec(chull,points,threads=3):
    pool = CHullCheck_pool(chull,threads=threads) 
    out = pool.map(CHullCheck, points)
    pool.close()
    pool.join()
    del pool
    return out

if __name__ == "__main__":


    ########################

    topology = aqtests.get("simple.prmtop")
    trajectory = aqtests.get("simple.nc")

    ########################
    
    log.message("Read trajectory...")
    reader = ReadAmberNetCDFviaMDA(topology, trajectory)

    traj_object = "(resname WAT) and (around 6 (resnum 88 or resnum 90 or resnum 136))"
    traj_scope = "protein"
    
    scope = reader.parse_selection(traj_scope)
    
    max_frame = float('inf')
    max_frame = 999

    ########################

    log.message("Loop over frames: search of waters in object...")
    pbar = log.pbar(max_frame)
    
    # loop over frames
    # scan for waters in object
    waters_ids_in_object_over_frames = {}
    all_H2O= None
    
    for frame in reader.iterate_over_frames():
        if frame > max_frame:
            break

        # current water selection
        H2O = reader.parse_selection(traj_object)
        
        # find convex hull of protein
        chull = scope.get_convexhull_of_atom_positions()

        H2O_coords = list(H2O.center_of_mass_of_residues())
        is_wat_within_chull = CHullCheck_exec(chull, H2O_coords)
        
        # discard wat out of scope
        H2O_new = None
        for iwwc,wat in zip(is_wat_within_chull,H2O.iterate_over_residues()):
            if iwwc:
                if H2O_new is None:
                    H2O_new = wat
                else:
                    H2O_new += wat
        
        # add it to all waters in object
        if all_H2O:
            all_H2O += H2O_new
            all_H2O.uniquify()
        else:
            all_H2O = H2O_new
        
        # remeber ids of water in object in current frame
        waters_ids_in_object_over_frames.update({frame:H2O_new.unique_resids()})
        
        pbar.update(frame)

    pbar.finish()
        
    log.message("Number of residues to trace: %d" % all_H2O.unique_resids_number())
        
    log.message("Init paths container...")
    # type and frames, consecutive elements correspond to residues in all_H2O
    paths = [GenericPath() for dummy in xrange(all_H2O.unique_resids_number())]
  
    log.message("Trajectory scan...")
    pbar = log.pbar(max_frame)

    # loop over frames
    for frame in reader.iterate_over_frames():
        if frame > max_frame:
            break
        
        # find convex hull of protein
        chull = scope.get_convexhull_of_atom_positions()

        all_H2O_coords = list(all_H2O.center_of_mass_of_residues())
        is_wat_within_chull = CHullCheck_exec(chull, all_H2O_coords)

        all_resids = [wat.first_resid() for wat in all_H2O.iterate_over_residues()]
        

        for nr,cir in enumerate(zip(all_H2O_coords,is_wat_within_chull,all_resids)):
            coord,ischull,resid = cir
            paths[nr].add_coord(coord)
            if ischull:
                # in scope
                if resid not in waters_ids_in_object_over_frames[frame]:
                    paths[nr].add_scope(frame)
                else:
                    # in object
                    paths[nr].add_object(frame)

        pbar.update(frame)    
    
    pbar.finish()
    
    log.message("Discard residues with empty paths...")
    
    
    zipped = [(wat,path) for wat,path in zip(all_H2O.iterate_over_residues(),paths) if len(path.frames) > 0]
    all_H2O,paths = zip(*[(wat,path) for wat,path in zip(all_H2O.iterate_over_residues(),paths) if len(path.frames) > 0])
    
    
            
    '''
    print "Extract coordinates..."
    
    for nr,wat in enumerate(all_H2O.iterate_over_residues()):
        coords = [wat.center_of_mass() for frame in reader.iterate_over_frames() if frame <= max_frame]
        print coords
    ''' 
    
    
    exit(0)


    '''
    import itertools

    # http://stackoverflow.com/questions/4628333/converting-a-list-of-integers-into-range-in-python
    def ranges(i):
        for a, b in itertools.groupby(enumerate(i), lambda (x, y): y - x):
            b = list(b)
            lo,hi = b[0][1], b[-1][1]
            if lo == hi:
                yield "resnum %d" % lo
            else:
                yield "resnum %d:%d" % (lo,hi)

    traj_it = " or ".join(ranges(all_waters))
    '''

    print "create all_H2O selection"
    all_H2O.uniquify()
    print all_H2O
    print waters_ids_in_object_over_frames

    for frame in reader.iterate_over_frames():
        if frame > max_frame:
            break
        print frame
        # get chull
        chull = scope.get_convexhull_of_atom_positions()
        print "convexhull points no %d" % len(chull.vertices_points)
        for wat in all_H2O.iterate_over_residues():
            if wat.resids.tolist()[0] not in waters_ids_in_object_over_frames[frame].tolist():
                if wat.resids.tolist()[0] in waters_ids_in_scope_over_frames[frame].tolist():
                    if chull.point_within(wat.center_of_mass()):
                        pass
                    #print "%r in scope" % wat,
        #print ""

