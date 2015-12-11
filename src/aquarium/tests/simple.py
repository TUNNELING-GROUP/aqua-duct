
from aquarium import tests as aqtests
from aquarium.traj.reader import ReadAmberNetCDFviaMDA

if __name__ == "__main__":

    topology = aqtests.get("simple.prmtop")
    trajectory = aqtests.get("simple.nc")

    reader = ReadAmberNetCDFviaMDA(topology, trajectory)

    # TEST AREA -->

    print reader.topology_file_name
    print reader.trajectory_file_name
    print reader.trajectory_object
    print reader.number_of_frames
    #print list(reader.iterate_over_frames())
    print "CA"
    CA = reader.parse_selection("protein and name CA")
    for frame in reader.iterate_over_frames():
        if frame >= 5:
            break
        print frame,CA.center_of_mass()

    print "H2O"
    H2O = reader.parse_selection("(resname WAT) and (around 6 (resnum 88 or resnum 90 or resnum 136))")
    for frame in reader.iterate_over_frames():
        H2O = reader.parse_selection("(resname WAT) and (around 6 (resnum 88 or resnum 90 or resnum 136))")
        if frame >= 5:
            break
        print H2O.unique_resids()
        print list(H2O.center_of_mass_of_residues())
    
    print "convexhull"
    print H2O.get_convexhull_of_atom_positions()
    print '\n'.join(dir(H2O.get_convexhull_of_atom_positions()))
    
    # <-- TEST AREA
    
    print "real data"
    
    traj_scope = "protein"

    max_frame = float('inf')
    #max_frame = 19
    traj_scope = "protein"
    traj_over_scope = "(resname WAT) and (around 2 protein)"
    traj_object = "(resname WAT) and (around 6 (resnum 88 or resnum 90 or resnum 136))"

    scope = reader.parse_selection(traj_scope)
    # find all waters that are in object at any point in the simulation
    waters_ids_in_object_over_frames = {}
    waters_ids_in_scope_over_frames = {}

    all_H2O = None
    for frame in reader.iterate_over_frames():
        if frame > max_frame:
            break
        print frame
        H2O = reader.parse_selection(traj_object)
        H2O_in_scope = reader.parse_selection(traj_over_scope)
        if all_H2O:
            all_H2O += H2O
            all_H2O.uniquify()
        else:
            all_H2O = H2O
        print all_H2O
        waters_ids_in_object_over_frames.update({frame:H2O.unique_resids()})
        waters_ids_in_scope_over_frames.update({frame:H2O_in_scope.unique_resids()})


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

