#!/usr/bin/env python

#from distutils.core import setup
from setuptools import find_packages,setup

version = __import__('aqueduct').version() # version tuple
version = '.'.join(map(str,version)) # version str

setup(name='aqueduct',
      version=version,
      description='Tracing residues in MD simulation',
      author='Tomasz Magdziarz',
      author_email='tomasz.magdziarz@polsl.pl',
      url='http://tunnelinggroup.pl',
      packages=find_packages(include=['aqueduct*']),
      install_requires=['numpy>=1.7',
                        'scipy>=0.13',
                        'scikit-learn>=0.14',
                        'MDAnalysis[AMBER]>=0.14',
                       ],
      extras_requires={'visual': ['matplotlib','pymol>=1.4'],
                       'valve': ['roman']},
      scripts=['apps/valve'],
      provides=['aqueduct'],
      classifiers=['Development Status :: 4 - Beta',
                   'Environment :: Console',
                   'Intended Audience :: Bioinformaticians',
                   'Operating System :: POSIX',
                   'Programming Language :: Python :: 2.7',
                  ],
      )
