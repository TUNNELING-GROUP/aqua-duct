#!/bin/sh

# install necessary libraries - required by netCDF4
apt-get install libnetcdf-dev libhdf5-dev

PIP=`which pip`

# clear LD_LIBRARY_PATH <-- customize it to your system requirements
export LD_LIBRARY_PATH=

$PIP install "netCDF4>=1.0" $@

$PIP install "MDAnalysis[AMBER]>=0.14" $@
