#!/usr/bin/env python2
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

import sys
from setuptools import find_packages, setup

version = __import__('aquaduct').version_nice()  # version str

install_requires =['numpy>=1.7.0',
                  'scipy>=0.14.0',
                  'scikit-learn>=0.16.0',
                  'MDAnalysis[amber]==0.15.0',
                  'roman>=2.0.0',
                  ]

def install_requires_nice(level=0):
      import re
      for ir in install_requires:
            print ("    "*level) + "* " + " ".join(re.split('(>=|==|<=|>|<|=)',ir))

author = __import__('aquaduct').__author__  # version str


setup(name='aquaduct',
      version=version,
      description='Tracing residues in MD simulation',
      author=author,
      author_email='info@aquaduct.pl',
      url='http://aquaduct.pl',
      packages=find_packages(include=['aquaduct*']),
      install_requires=install_requires,
      extras_require={'full_pymol': ["pymol>=1.4"],
                      'graphs': ['matplotlib'],
                      },
      scripts=['apps/valve.py', 'apps/valve_run'],
      provides=['aquaduct'],
      classifiers=['Development Status :: 4 - Beta',
                   'Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'Operating System :: POSIX',
                   'Programming Language :: Python :: 2.7',
                   ],
      )
