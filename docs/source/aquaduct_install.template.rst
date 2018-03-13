Aqua-Duct installation guide
============================

Overview
--------

Aqua-Duct software is software written in Python (CPython) and comprises of two elements:

#. aquaduct - a Python package,
#. valve - a script that uses :mod:`aquaduct` to perform calculations.

**Download**

You can download Aqua-Duct packages directly from `Aqua-Duct homepage <http://aquaduct.pl>`_. This page includes older versions of Aqua-Duct as well as development version.

If you follow this installation guide you will install current release.


Troubleshooting
---------------

If you encounter any problems with installation do not hesitate to contact us at `info@aquaduct.pl <info@aquaduct.pl>`_. We are **REALLY** willing to help!

Please, provide us with us much info as you can. In particular try to include following information:

* Operating system's name and version, and CPU architecture (if relevant).
* Python version.
* Command(s) you have used for installation.
* Any error/warning/info message(s) that emerged during or after installation.

Requirements
------------

Software-wise requirements
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. include:: aquaduct_install_requires.rst

Hardware-wise requirements
^^^^^^^^^^^^^^^^^^^^^^^^^^

Aqua-Duct should work on every machine on which you can install the above mentioned software. On computers older than 10 years it may work very slow though. We recommend 64bit SMP architecture, with at least 4GB RAM (32 GB RAM is recommended).

Installation
------------

Generic Python installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The easiest way to install Aqua-Duct is to install Python 2.7 and use following command::

    pip AQPIP

If *pip* is not available try to install it by typing::

    easy_install pip

Depending on the settings of your system you can prepend the above command with `sudo` or `doas` or do *user* installation::

    # sudo
    sudo pip AQPIP

    # doas
    doas pip AQPIP

    # 'user' installation
    pip AQPIP --user

It is also good idea to try to install Aqua-Duct using virtualenv::

    virtualenv aquaduct_installation
    cd aquaduct_installation
    . bin/activate
    pip AQPIP

Installation of PyMOL
#####################

Under most modern GNU/Linux distributions PyMOL is available as a package in repositories. For example if you are under Ubuntu/Debian you can install it by following command::

    sudo apt-get install pymol

Under Windows there are several ways to install PyMOL, for more details see `PyMOL web site <http://pymol.org>`_.

Instructions for macOS and OpenBSD are in appropriate sections below.

GNU/Linux
^^^^^^^^^

Installation was tested on limited number of GNU/Linux systems. On the most of modern installations you can simply follow generic instructions, for example under Ubuntu 16.04 you can type::

    sudo pip AQPIP

NetCDF4 & MDAnalysis installation Ubuntu 14.04
##############################################

Other systems may require additional work, in particular installation of NetCDF4 is sometimes cumbersome. Following is an example how to install all required packages under Ubuntu 14.04::

    # install required python packages
    sudo apt-get install python-dev python-pip python-numpy python-scipy python-matplotlib python-scikits-learn

    # install necessary libraries and git -  all required to compile netCDF4
    sudo apt-get -y install libnetcdf-dev libhdf5-dev git

    # clone netcdf4 python repository
    git clone https://github.com/Unidata/netcdf4-python.git
    # cd to cloned repository
    cd netcdf4-python
    # modify setup.cfg to add paths of hdf5 and netcdf4 libraries
    sed -i '/\[directories\]/a \
    HDF5_dir = /usr/lib \
    HDF5_libdir = /usr/lib \
    HDF5_incdir = /usr/include \
    netCDF4_dir = /usr/lib \
    netCDF4_libdir = /usr/lib \
    netCDF4_incdir = /usr/include' setup.cfg
    # run setup.py
    sudo python setup.py install

    # install MDAnalysis
    sudo pip install "MDAnalysis[amber]==0.16.2"

If everything went fine you can follow generic instructions.

SciPy update and Ubuntu/Debian
##############################

Debian (and Ubuntu) uses strange approach to Python installation. To install newer version of SciPy (if required) try following procedure::

    # install libraries required for SciPy compilation
    apt-get build-dep python-scipy

    # install SciPy
    easy_install-2.7 --upgrade scipy

.. warning::

    The above procedure will remove current SciPy from `easy-install.pth` file.

macOS
^^^^^

Aqua-Duct installation was tested on macOS Sierra and is quite straightforward. It can be installed either with existing system Python or with custom Python installation. In both cases one have to install Xcode for the App Store.

System native Python
####################

::

    sudo easy_install pip
    sudo pip AQPIP

The drawback of using system Python installation is a lack of PyMOL. It should be, however, relatively easy to compile PyMOL on your own. Try to follow compilation instruction under BSD systems.

Custom Python
#############

This is recommended way of Aqua-Duct installation. If you do not have custom Python installation you can get it by using one of package managers available for macOS, for example `homebrew <http://brew.sh/>`_. With this package manager you can do following::

    brew install python
    sudo easy_install pip
    sudo pip AQPIP

Next, you can install PyMOL::

    brew install pymol
    brew cask install xquartz

Once XQuartz is installed you should reboot. The above procedure installs PyMOL, however, PyMOL Python modules are not visible. To fix it you can issue following commands::

    cd /usr/local/lib/python2.7/site-packages
    sudo ln -s /usr/local/Cellar/pymol/*/libexec/lib/python2.7/site-packages/* ./

The above instruction assumes that you are using `brew` and you have only one PyMOL installation.

Windows
^^^^^^^

Installation under Windows is also possible. The limiting factor is MDAnalysis which is not officially available under Windows yet. You can, however, install Cygwin and perform Aqua-Duct installation in Cygwin.

First, start with `Cygwin installation <https://cygwin.com/>`_. During the setup select following packages:

* python (2.7)
* python-devel (2.7)
* python-cython
* libnetcdf-devel
* libhdf5-devel
* liblapack-devel
* libopenblas
* python-numpy
* python-six

Another key component that have to be installed is C, C++ and Fortran compilers. You can simply install **gcc-g++** and **gcc-fortran** packages as a first choice, select following packages:

* gcc-g++
* gcc-fortran

Once Cygwin is installed with all required libraries you can perform following steps::

    # install pip
    easy_install-2.7 pip

First, try to install SciPy::

    # install SciPy
    pip install scipy

If you encounter any problems related to missing **xlocale.h** header file try the following workaround::

    # prepare fake xlocale.h
    ln -s /usr/include/locale.h xlocale.h
    export CFLAGS="I"$( pwd )

    # install SciPy
    pip install scipy

.. note::

	The above procedure for SciPy installation might not be optimal. For more information please got to `SciPy web page <https://www.scipy.org/>`_.

Now, install **scikit-learn** and then Aqua-Duct::

    # install scikit-learn
    pip install scikit-learn

    # finally, install aquaduct
    pip AQPIP

OpenBSD
^^^^^^^

Aqua-Duct can be also installed under OpenBSD (5.9 and 6.0 amd64).
NetCDF-c version 4 has to be installed as OpenBSD ships only netCDF in version 3. First, install hdf5 library and GNU make::

    # install hdf5 and GNU make
    pkg_add hdf5 gmake

Next, download netCDF sources. Version 4.2.1.1 works out of the box but is a bit outdated. Visit `NetCDF web page <https://www.unidata.ucar.edu/software/netcdf/>`_ and select version of your choice. Older versions are available in the `FTP archive <ftp://ftp.unidata.ucar.edu/pub/netcdf/old/>`_. Once netCDF is downloaded and extracted go to the source directory and try following procedure::

    # set LD and CPP flags
    export LDFLAGS=-L/usr/local/lib
    export CPPFLAGS=-I/usr/local/include

    # configure project
    ./configure --enable-shared --enable-dap --disable-doxygen --enable-netcdf-4 --prefix=/path/to/netCDF4/lib

    # make and install
    gmake
    gmake install

You may now install py-scipy package::

    pkg_add py-scipy

Install pip if it is missing::

    pkg_add py-pip

Install netCDF4 Python::

    # define netcdf-4 installation directory
    export NETCDF4_DIR=/path/to/netCDF4/lib
    pip2.7 install netCDF4

At this point you can follow generic Python instructions, type::

    pip2.7 AQPIP

PyMOL at OpenBSD
################

According to our knowledge it is possible to install PyMOL 1.4.1 and it is sufficient to work with Aqua-Duct. Go to `SourceForge PyMOL download page <https://sourceforge.net/projects/pymol/files/pymol/>`_ and download, save, and extract sources.

PyMOL requires Python Mega Widgets. Download, for example Pmw 1.3.3b from `SourceForge Pmw download page <https://sourceforge.net/projects/pmw/files/Pmw/>`_. Extract it and install by::

    python2.7 setup.py install

TKinter (2.7) and several other packages are also required::

    pkg_add python-tkinter freeglut glew png

Next, go to the extracted PyMOL sources open setup.py and modify inc_dirs variable at line 129 by adding following paths::

    "/usr/X11R6/include/freetype2",
    "/usr/X11R6/include",
    "/usr/local/include",

Now, you can build and install PyMOL by typing following commands::

    python2.7 setup.py build
    python2.7 setup.py install
    python2.7 setup2.py install
    cp pymol /usr/local/bin

PyMOL can be run by typing `pymol` or can be used as Python module.

Other BSDs
##########

Installation on other BSDs might be easier. For example, Python netCDF4 is available in ports of FreeBSD and DragonFlyBSD. Try to install it and SciPy, then proceed to generic Python installation instructions.

If you are using NetBSD or other BSD try to follow OpenBSD instructions.

