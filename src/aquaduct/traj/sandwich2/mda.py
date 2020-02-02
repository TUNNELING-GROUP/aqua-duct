# -*- coding: utf-8 -*-



from aquaduct.traj.sandwich2.reader import BaseReader

from aquaduct import logger
from aquaduct.utils.helpers import version_parser
from os.path import splitext
import re

import MDAnalysis as mda

################################################################################
# Check MDAnalysis version

def mda_ver():
    return version_parser(mda.__version__)


# FIXME: do it according to user options
if mda_ver() < version_parser('0.16'):
    logger.error('Unsupported MDAnalysis version %s; should be 0.16.2 or > 0.19.', mda.__version__)
    raise NotImplementedError('Unsupported MDAnalysis version %s; should be 0.16.2 or > 0.19.' % mda.__version__)

if mda_ver() < version_parser('0.17') and mda_ver() < version_parser('0.20'):
    logger.warning('Unsupported MDAnalysis version %s.', mda.__version__)
    logger.warning('Trying to mitigate potential problems by setting `use_periodic_selections = False`.')
    mda.core.flags['use_periodic_selections'] = False


################################################################################
# raw

mda_available_formats = {re.compile('(nc|NC)'): 'nc',
                         re.compile('(prmtop|parmtop|top|PRMTOP|PARMTOP|TOP)'): 'PRMTOP',
                         re.compile('(dcd|DCD)'): 'LAMMPS',
                         re.compile('(psf|PSF)'): 'psf',
                         re.compile('(pdb|PDB)'): 'pdb',
                         re.compile('(crd|CRD)'): 'crd',
                         re.compile('(xtc|XTC)'): 'XTC'}

def open_raw(topology, trajectory):
    topology_ext = splitext(topology)[1][1:]
    for afk in list(mda_available_formats.keys()):
        if afk.match(topology_ext):
            topology_ext = mda_available_formats[afk]
            break
    trajectory_ext = splitext(trajectory[0])[1][1:]
    for afk in list(mda_available_formats.keys()):
        if afk.match(trajectory_ext):
            trajectory_ext = mda_available_formats[afk]
            break
    return mda.Universe(topology,
                        trajectory,
                        topology_format=topology_ext,
                        format=trajectory_ext)


################################################################################


class Reader(BaseReader):

    def open_trajectory(self):
        # returns raw trajectory objet to be interpreted by this class
        return open_raw(self.topology,self.trajectory)


    def close_trajectory(self):
        if hasattr(self, 'trajectory_object'):
            if hasattr(self.trajectory_object, 'trajectory'):
                if hasattr(self.trajectory_object.trajectory, 'close'):
                    self.trajectory_object.trajectory.close()

    def physical_number_of_frames(self):
        return len(self.trajectory_object.trajectory)
