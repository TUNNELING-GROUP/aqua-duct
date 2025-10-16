#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2025
#  

from setuptools import setup


setup(
    name = "aquaduct",
    python_requires = ">=3.10",
    package_dir = {"": "src"},
    keywords = "molecular-dynamics solvent",
    scripts = [
        "src/aquaduct/scripts/valve.py",
        "src/aquaduct/scripts/valve_run",
        "src/aquaduct/scripts/valveconfig.py",
        "src/aquaduct/scripts/valveconfig_run",
        "src/aquaduct/scripts/portal.py",
        "src/aquaduct/scripts/portal_run",
        "src/aquaduct/scripts/pond.py",
        "src/aquaduct/scripts/pond_run",
        "src/aquaduct/scripts/kraken.py",
        "src/aquaduct/scripts/kraken_run",
        "src/aquaduct/scripts/hs_resize.py"
     ],
    classifiers = [
        "Development Status :: 3 - Alpha",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Operating System :: POSIX",
        "Programming Language :: Python"
        ],
    package_data = {"": ["apps/valveconfig/*.gif"]},
    include_pakacge_data = True
    )
