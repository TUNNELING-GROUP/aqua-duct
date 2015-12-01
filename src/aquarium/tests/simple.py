
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
    CA = reader.parse_selection("protein and name CA")
    for frame in reader.iterate_over_frames():
        if frame >= 5:
            break
        print frame,CA.center_of_mass()

    H2O = reader.parse_selection("(resname WAT) and (around 6 (resnum 88 or resnum 90 or resnum 136))")

    for frame in reader.iterate_over_frames():
        H2O = reader.parse_selection("(resname WAT) and (around 6 (resnum 88 or resnum 90 or resnum 136))")
        if frame >= 5:
            break
        print H2O.unique_resnums()
        print list(H2O.center_of_mass_of_residues())

    print H2O.get_convexhull_of_atom_positions()
    print '\n'.join(dir(H2O.get_convexhull_of_atom_positions()))
    # <-- TEST AREA

    traj_scope = "protein"
    traj_object = "(resname WAT) and (around 6 (resnum 88 or resnum 90 or resnum 136))"

    # find all waters that are in object at any point in the simulation
    all_waters = []
    for frame in reader.iterate_over_frames():
        print frame
        H2O = reader.parse_selection(traj_object)
        all_waters.extend(H2O.unique_resnums())
    print set(all_waters)
