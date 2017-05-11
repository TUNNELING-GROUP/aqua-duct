# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2016-2017  Tomasz Magdziarz, Alicja Płuciennik, Michał Stolarczyk <info@aquaduct.pl>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import logging
logger = logging.getLogger(__name__)
from os import unlink

import MDAnalysis as mda

from aquaduct.utils import clui
from aquaduct.utils.helpers import create_tmpfile


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
                reader.set_real_frame(frame)
                self.mdawriter.write(to_dump)
            else:
                logger.error(
                    'Requested frame %d exceeded available number of frames %d.' % (frame, reader.number_of_frames - 1))

    def close(self):
        self.mdawriter.close()
        return self.pdbfile

    def __del__(self):
        unlink(self.pdbfile)
