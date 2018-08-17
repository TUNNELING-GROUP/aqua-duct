# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2016-2018  Tomasz Magdziarz, Alicja Płuciennik, Michał Stolarczyk <info@aquaduct.pl>
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

from aquaduct import logger
import os

class WritePDB(object):

    def __init__(self,pdbfile,csvfile=None,scale_bf=0.1):

        self.current_atom = 1
        self.scale_bf = scale_bf

        self.fh = open(pdbfile,'w')
        self.fh_csv = None
        if csvfile:
            self.fh_csv = open(csvfile,'w')

    def print_atom_line(self,xyz,bf):
        atom = "%5d" % self.current_atom
        x = ((" "*8)+("%0.3f" % xyz[0])[:8])[-8:]
        y = ((" "*8)+("%0.3f" % xyz[1])[:8])[-8:]
        z = ((" "*8)+("%0.3f" % xyz[2])[:8])[-8:]
        #x = ("        %0.3f" % xyz[0])[-8:]
        #y = ("        %0.3f" % xyz[1])[-8:]
        #z = ("        %0.3f" % xyz[2])[-8:]
        b = ((" "*4)+("%0.2f" % (bf*self.scale_bf))[:4])[-4:]
        return ("ATOM  %s  H   FIL T   1    %s%s%s        %s" % (atom,x,y,z,b))+os.linesep

    def print_conect_line(self,a1,a2):
        atom1 = "%5d" % a1
        atom2 = "%5d" % a2
        return ("CONECT%s%s" % (atom1,atom2))+os.linesep

    def write_connected(self, line, bf):
        conect = None
        for xyz,b in zip(line, bf):
            self.fh.write(self.print_atom_line(xyz,b))
            if conect:
                self.fh.write(self.print_conect_line(conect,self.current_atom))
            elif self.fh_csv is not None:
                self.fh_csv.write("X,Y,Z"+os.linesep) # header
            if self.fh_csv is not None:
                self.fh_csv.write(("%f,%f,%f"+os.linesep) % tuple(xyz))
            conect = self.current_atom
            self.current_atom += 1

    def write_scatter(self, scatter, bf):
        for xyz,b in zip(scatter, bf):
            self.fh.write(self.print_atom_line(xyz,b))
            self.current_atom += 1

    def __enter__(self):
        return self

    def __exit__(self, typ, value, traceback):
        if typ is None:
            self.__del__()


    def __del__(self):
        self.fh.close()
        if self.fh_csv is not None:
            self.fh_csv.close()
