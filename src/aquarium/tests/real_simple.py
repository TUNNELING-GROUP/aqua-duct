
from aquarium import tests as aqtests
from aquarium.traj.reader import ReadAmberNetCDFviaMDA
from aquarium.traj.paths import GenericPath

if __name__ == "__main__":

    topology = aqtests.get("simple.prmtop")
    trajectory = aqtests.get("simple.nc")

    print "Read trajectory..."
    reader = ReadAmberNetCDFviaMDA(topology, trajectory)

    traj_object = "(resname WAT) and (around 6 (resnum 88 or resnum 90 or resnum 136))"
    traj_scope = "protein"
    
    scope = reader.parse_selection(traj_scope)
    
    max_frame = float('inf')
    max_frame = 19

    print "Loop over frames: search of waters in object..."
    # loop over frames
    # scan for waters in object
    waters_ids_in_object_over_frames = {}
    all_H2O= None
    for frame in reader.iterate_over_frames():
        if frame > max_frame:
            break
        print "Frame",frame

        # current water selection
        H2O = reader.parse_selection(traj_object)
        # add it to all waters in object
        if all_H2O:
            all_H2O += H2O
            all_H2O.uniquify()
        else:
            all_H2O = H2O
        # remeber ids of water in object in current frame
        waters_ids_in_object_over_frames.update({frame:H2O.unique_resids()})
        
        
    print "Init paths container..."
    # type and frames, consecutive elements correspond to residues in all_H2O
    paths = [GenericPath() for dummy in xrange(all_H2O.unique_resids_number())]

    print "Trajectory scan..."
    # loop over frames
    for frame in reader.iterate_over_frames():
        if frame > max_frame:
            break
        print "Frame",frame
        
        # find convex hull of protein
        chull = scope.get_convexhull_of_atom_positions()
        print "\tCurrent %s convexhull of %d points" % (traj_scope, len(chull.vertices_ids))
        # loop over waters in all waters
        for nr,wat in enumerate(all_H2O.iterate_over_residues()):
            if chull.point_within(wat.center_of_mass()):
                # in scope
                if wat.first_resid() not in waters_ids_in_object_over_frames[frame]:
                    paths[nr].add_scope(frame)
                else:
                    # in object
                    paths[nr].add_object(frame)
            
    print "Extract coordinates..."
    
    for nr,wat in enumerate(all_H2O.iterate_over_residues()):
        coords = [wat.center_of_mass() for frame in reader.iterate_over_frames() if frame <= max_frame]
        print coords
        
    
    
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

