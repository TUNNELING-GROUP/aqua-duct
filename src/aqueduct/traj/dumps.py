# -*- coding: utf-8 -*-


from os import unlink

import MDAnalysis as mda

from aqueduct.utils import log
from aqueduct.utils.helpers import create_tmpfile


class TmpDumpWriterOfMDA(object):
    def __init__(self):

        # creates tmp file
        self.pdbfile = create_tmpfile(ext='pdb')

        # create mda writer
        self.mdawriter = mda.Writer(self.pdbfile, multiframe=True)

    def dump_frames(self, reader, frames, selection='protein'):

        to_dump = reader.parse_selection(selection)

        for frame in frames:
            if frame < reader.number_of_frames:
                reader.set_current_frame(frame)
                self.mdawriter.write(to_dump)
            else:
                log.error(
                    'Requested frame %d exceeded available number of frames %d.' % (frame, reader.number_of_frames - 1))

    def close(self):
        self.mdawriter.close()
        return self.pdbfile

    def __del__(self):
        unlink(self.pdbfile)
