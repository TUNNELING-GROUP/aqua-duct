Aqua-Duct installation guide
===========================

Overview
--------


This package comprises of two elements:

#. aquaduct,
#. valve.

Aqua-Duct is a Python module. It is a collection of tools to trace residues in MD simulation. Valve is a driver Python script. It uses aquaduct to perform such a tracing. Valve is shipped together with Aqua-Duct.

Install
-------

Aqua-Duct
^^^^^^^^

Installation was tested on limited number of POSIX-like systems.

In some specific cases installation is very simple:

#. Download :download:`aquaduct_0.2.25.tar.gz` bundle file.
#. Unpack aquaduct bundle file.
#. Go to src directory.
#. Type::

    python setup.py install

Aqua-Duct requires several Python modules to work and in particular it requires MDAnalysis with AMBER support. This, on the other hand, requires netCDF4. Installation of this combination is sometimes cumbersome. General procedure is following:

#. Install libnetcdf4 and libhdf5 development libraries.
#. Install netCDF4::

    pip install netCDF4

#. Try to install aquaduct.

If, by chance, you are on Ubuntu 14.04 you can try helper script :download:`ubuntu_mdanalysis_install_helper.sh`.

Valve
^^^^^

Valve does not need installation *per se*. Once aquaduct is installed, valve can be run by a following command::

    valve.py --help

Valve script, ie valve.py, is located in apps directory.

Extras
^^^^^^

Access to some visualization capabilities of Aqua-Duct requires additional Python modules:

#. matplotlib,
#. pymol.

These are usually easy to install.

Troubleshooting
---------------

If you encounter any problems with installation do not hesitate to contact us at `<info@aquaduct.pl>`_. We are REALLY willing to help!
