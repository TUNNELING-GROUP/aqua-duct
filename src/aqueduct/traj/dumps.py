

import MDAnalysis as mda
from aqueduct.utils.helpers import create_tmpfile
from os import unlink

class TmpDumpWriterOfMDA(object):


    def __init__(self):

        # creates tmp file
        self.pdbfile = create_tmpfile(ext='pdb')

        # create mda writer
        self.mdawriter = mda.Writer(self.pdbfile, multiframe=True)

    def dump_frames(self,reader,frames,selection='protein'):

        to_dump = reader.parse_selection(selection)

        for frame in frames:
            reader.set_current_frame(frame)
            self.mdawriter.write(to_dump)

    def close(self):
        self.mdawriter.close()
        return self.pdbfile

    def __del__(self):
        unlink(self.pdbfile)


