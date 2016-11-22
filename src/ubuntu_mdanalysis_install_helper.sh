#!/bin/sh

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2016  Tomasz Magdziarz, Alicja Płuciennik, Michał Stolarczyk <info@aquaduct.pl>
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

# install necessary libraries - required by netCDF4
sudo apt-get -y install libnetcdf-dev libhdf5-dev git

# install required python packages
sudo apt-get -y install python-dev python-pip python-numpy python-scipy python-matplotlib python-scikits-learn pymol

PIP=`which pip`

# clear LD_LIBRARY_PATH <-- customize it to your system requirements
export LD_LIBRARY_PATH=

# do git install of netCDF4
BUILDDIR=`mktemp -d`
CWD=`pwd`
cd $BUILDDIR
git clone https://github.com/Unidata/netcdf4-python.git
cd netcdf4-python
sed -i '/\[directories\]/a \
HDF5_dir = /usr/lib \
HDF5_libdir = /usr/lib \
HDF5_incdir = /usr/include \
netCDF4_dir = /usr/lib \
netCDF4_libdir = /usr/lib \
netCDF4_incdir = /usr/include' setup.cfg
sudo python setup.py install
cd $CWD
rm -rf $BUILDDIR

# install MDAnalysis
sudo $PIP install "MDAnalysis[amber]>=0.15" $@
