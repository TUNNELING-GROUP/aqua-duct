#!/bin/sh

# install necessary libraries - required by netCDF4
apt-get -y install libnetcdf-dev libhdf5-dev

# install required python packages
apt-get -y  python-pip python-numpy python-scipy python-matplotlib python-scikit-learn pymol

PIP=`which pip`

# clear LD_LIBRARY_PATH <-- customize it to your system requirements
export LD_LIBRARY_PATH=

# install netCDF4 and MDAnalysis
$PIP install "netCDF4>=1.0" $@
$PIP install "MDAnalysis[amber]>=0.15" $@
