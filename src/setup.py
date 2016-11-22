#!/usr/bin/env python2
# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2016  Tomasz Magdziarz <info@aquaduct.pl>
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

import sys
from setuptools import find_packages, setup

version = __import__('aqueduct').version()  # version tuple
version = '.'.join(map(str, version))  # version str
#sys.path.append('apps')
#version = '.'.join((version, __import__('valve').version_onenumber()))

setup(name='aqueduct',
      version=version,
      description='Tracing residues in MD simulation',
      author='Tomasz Magdziarz',
      author_email='tomasz.magdziarz@polsl.pl',
      url='http://tunnelinggroup.pl',
      packages=find_packages(include=['aqueduct*']),
      install_requires=['numpy>=1.7',
                        'scipy>=0.13',
                        'scikit-learn>=0.16',
                        'MDAnalysis[amber]>=0.15',
                        'roman',
                        ],
      extras_require={'full_pymol': ["pymol>=1.4"],
                      'graphs': ['matplotlib'],
                      },
      scripts=['apps/valve.py', 'apps/valve_run'],
      provides=['aqueduct'],
      classifiers=['Development Status :: 4 - Beta',
                   'Environment :: Console',
                   'Intended Audience :: Bioinformaticians',
                   'Operating System :: POSIX',
                   'Programming Language :: Python :: 2.7',
                   ],
      )
