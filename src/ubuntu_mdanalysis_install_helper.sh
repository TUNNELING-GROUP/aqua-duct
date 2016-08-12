#!/bin/sh

# install necessary libraries - required by netCDF4
apt-get -y install libnetcdf-dev libhdf5-dev git

# install required python packages
apt-get -y install python-dev python-pip python-numpy python-scipy python-matplotlib python-scikits-learn pymol

PIP=`which pip`

# clear LD_LIBRARY_PATH <-- customize it to your system requirements
export LD_LIBRARY_PATH=

# do git install of netCDF4
BUILDDIR=`mktemp -d`
CWD=`pwd`
cd $BUILDDIR
git clone https://github.com/Unidata/netcdf4-python.git
cd netcdf4-python
sed '/\[directories\]/a \
HDF5_dir = /usr/lib \
HDF5_libdir = /usr/lib \
HDF5_incdir = /usr/include \
netCDF4_dir = /usr/lib \
netCDF4_libdir = /usr/lib \
netCDF4_incdir = /usr/include' setup.cfg
python setup.py install
cd $CWD
rm -rf $BUILDDIR

# install MDAnalysis
$PIP install "MDAnalysis[amber]>=0.15" $@
