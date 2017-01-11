Aqua-Duct installation guide
============================

Overview
--------

Aqua-Duct software is software written in Python (CPython) and comprises of two elements:

#. aquaduct - a Python package,
#. valve - a script that uses  :mod:`aquaduct` to perform calculations.

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

Aqua-Duct should work on every machine on which you can install the above mentioned software. On computers older then 10 years it may work very slow though. We recommend 64bit SMP architecture, with at least 4GB RAM (32 GB RAM is recommended).

Installation
------------

Generic Python installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The easiest way to install Aqua-Duct you have to install Python and use following command::

    pip install --extra-index-url https://testpypi.python.org/pypi aqueduct

Depending on the settings of your system you can prepend the above command with `sudo` or `doas` or do *user* instalation::

    # sudo
    sudo pip install --extra-index-url https://testpypi.python.org/pypi aqueduct

    # doas
    doas pip install --extra-index-url https://testpypi.python.org/pypi aqueduct

    # 'user' installation
    pip install --extra-index-url https://testpypi.python.org/pypi aqueduct --user

It is also good idea to try to install Aqua-Duct using virtualenv::

    virtualenv aquaduct_installation
    cd aquaduct_installation
    . bin/activate
    pip install --extra-index-url https://testpypi.python.org/pypi aqueduct

Installation of PyMOL
#####################

Under most modern GNU/Linux distributions PyMOL is available as a package in repositories. For example if you are under Ubuntu/Debian you can install it by following command::

    sudo apt-get install pymol

Under Windows there are several ways to install PyMOL, for more details see `PyMOL web site <http://pymol.org>`_.

Instructions for macOS are below.

GNU/Linux
^^^^^^^^^

Installation was tested on limited number of GNU/Linux systems. On the most of modern installations you can simply follow generic instructions, for example under Ubuntu 16.04 you can type::

    sudo pip install --extra-index-url https://testpypi.python.org/pypi aqueduct

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
    sudo pip install "MDAnalysis[amber]>=0.15"

If everything went fine you can follow generic instructions, type::

    sudo pip install --extra-index-url https://testpypi.python.org/pypi aqueduct

MacOS
^^^^^

Aqua-Duct installation was tested on MacOS Sierra and is quite straightforward. It can be installed either with existing system Python or with custom Python installation. In both cases one have to install Xcode for the App Store.

System native Python
####################

::

    sudo easy_install pip
    sudo pip install --extra-index-url https://testpypi.python.org/pypi aqueduct

The drawback of using system Python installation is a lack of PyMOL. It should be, however, realatively easy to compile PyMOL on your own. Try to follow compilation instruction under BSD systems.

Custom Python
#############

This is recomended way of aquaduct installation. If you do not have custom Python instalaltion you can get it by using one of package managers available for macOS, for example `homebrew <http://brew.sh/>`_. With this package manager you can do following::

    brew install python
    sudo easy_install pip
    sudo pip install --extra-index-url https://testpypi.python.org/pypi aqueduct

Next, you can install PyMOL::

    brew install pymol
    brew cask install xquartz

Once XQuartz is installed you should reboot. The above procedure installs PyMOL, however, PyMOL Python modules are not visible. To fix it you can issue following commands::

    cd /usr/local/lib/python2.7/site-packages
    sudo ln -s /usr/local/Cellar/pymol/*/libexec/lib/python2.7/site-packages/* ./

The above instruction assumes that you are using `brew` and you have only one PyMOL installation.

Windows
^^^^^^^

Installation under Windows is also possible. The limiting factor is MDAnalysis which is not officially available under Windows yet. You can, however, install Cygwin and perofrom Aqua-Duct installation in Cygwin.

First, start with `Cygwin installation <https://cygwin.com/>`_. During the setup select following packages:

* python (2.7)
* python-devel (2.7)
* python-cython
* libnetcdf-devel
* libhdf5-devel
* liblapack-devel
* libopenblas

You can also select following packages:

* python-numpy
* python-six

.. note::

    You can skip installation of these packages. If they are missing they will be installed automatically.

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
    pip install --extra-index-url https://testpypi.python.org/pypi aqueduct
