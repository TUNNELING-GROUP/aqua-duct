#!/usr/bin/env python2
# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2016-2018  Tomasz Magdziarz, Alicja Płuciennik, Michał Stolarczyk <info@aquaduct.pl>
# Copyright (C) 2018-2019  Tomasz Magdziarz, Michał Banas <info@aquaduct.pl>
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

from setuptools import find_packages, setup

version = __import__('aquaduct').version_nice()  # version str

install_requires = ['numpy>=1.10.0',  # this is required by MDA
                    'scipy>=0.17.1',
                    'scikit-learn>=0.16.0',
                    'MDAnalysis[amber]>=0.16.2',
                    'joblib>=0.13'
                    ]


def install_requires_nice(level=0):
    import re
    for ir in install_requires:
        print((" " * 4 * level) + "* " + " ".join(re.split('(>=|==|<=|>|<|=)', ir)))


author = __import__('aquaduct').__author__  # version str

with open('README.md', 'r') as fh:
    long_description = fh.read()

setup(name='aquaduct',
      version=version,
      description='Tracing molecules in MD simulation',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='http://aquaduct.pl',
      project_urls={
          'Documentation': 'https://tunneling-group.github.io/aqua-duct/',
          'Source': 'https://github.com/TUNNELING-GROUP/aqua-duct',
          'Tracker': 'https://github.com/TUNNELING-GROUP/aqua-duct/issues/',
      },
      author=author,
      author_email='info@aquaduct.pl',
      license='GNU GPL v3',
      keywords='molecular-dynamics solvent',
      packages=find_packages(include=['aquaduct*']),
      python_requires='>=3',
      install_requires=install_requires,
      extras_require={'full_pymol': ["pymol>=1.4"],
                      'graphs': ['matplotlib'],
                      },
      scripts=['apps/valve.py', 'apps/valve_run',
               'apps/valveconfig.py', 'apps/valveconfig_run',
               'apps/portal.py', 'apps/portal_run',
               'apps/pond.py', 'apps/pond_run',
               'apps/kraken.py', 'apps/kraken_run',
               'apps/hs_resize.py', ],
      provides=['aquaduct'],
      classifiers=['Development Status :: 5 - Production/Stable',
                   'Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'Operating System :: POSIX',
                   'Programming Language :: Python',
                   ],
      include_package_data=True
      )
