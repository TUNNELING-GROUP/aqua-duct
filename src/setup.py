#!/usr/bin/env python

from distutils.core import setup
from setuptools import find_packages



# Dynamically calculate the version
version = __import__('aqueduct').version_nice()

setup(name='Aqueduct',
      version=version,
      description='Tracing residues in MD simulation',
      author='Tomasz Magdziarz',
      author_email='tomasz.magdziarz@polsl.pl',
      url='http://tunnelinggroup.pl',
      packages=find_packages(include=['aqueduct*']),
      requires=['numpy (>=1.7)',
                'scipy (>=0.13)',
                'sklearn (>=0.14)',
                'MDAnalysis (>=0.12)',
                ],
      provides=['aqueduct'],
      classifiers=[
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'Intended Audience :: Bioinformaticians',
          'Operating System :: POSIX',
          'Programming Language :: Python :: 2.7',
          ],
      )
