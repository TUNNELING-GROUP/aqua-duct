
from aquarium import tests as aqtests
from aquarium.traj.reader import ReadAmberNetCDFviaMDA

if __name__ == "__main__":

    topology = aqtests.get("simple.prmtop")
    trajectory = aqtests.get("simple.nc")

    reader = ReadAmberNetCDFviaMDA(topology, trajectory)
    print reader.topology_file_name
    print reader.trajectory_file_name
    print reader.trajectory_object
    print reader.number_of_frames
    #print list(reader.iterate_over_frames())
    CA = reader.parse_selection("protein and name CA")
    for frame in reader.iterate_over_frames():
        print CA.center_of_mass()