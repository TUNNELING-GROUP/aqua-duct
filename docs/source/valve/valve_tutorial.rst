*Valve* tutorial
================

This is a tentative *Valve* manual. Created for the sake of Aqueduct training we have today. Eventually, it will be rewritten to the official version.

This tutorial assumes :mod:`aqueduct` and *Valve* is already installed - see :doc:`../aqueduct_install`. It is also assumed that user is acquientained with :doc:`valve_manual` and *Valve* :doc:`valve_config`.


*Valve* invocation
------------------

Usually *Valve* is run by::

    valve.py

Due to specific setup in our laboratory *Valve* has to be run through simple wrapper script::

    valve_run

Additionaly, to speed up all calculations it is assumed that *Valve* is run with ``--max-frame 1000`` option::

    valve_run --max-frame 1000

To check is *Valve* is isntalled and works properly try to issue following commands::

    valve_run --help
    valve_run --version

Test data
---------

**Mause!**

We will use 10ns Amber MD simulation data of sEH protein (PDBID **1cqz**). Necessary files can be downloaded `here <http://localhost:8001>`_:

* Go to download server.
* Go to ``1cqz`` directory.
* Download all files and save them in sane location on your machine. Please note, that ``.nc`` file is ca. 3.5 GB so it may take a while to download it.

Inspect your system
-------------------

Before we start any calculations lets have a look at the protein of interest. Start *PyMOL* and get ``1cqz`` PDB structure (for example by typing in *PyMOL* command prompt ``fetch 1cqa``).

To setup *Valve* calculations we need to know active site of the protein. More precisely we need to know IDs or residues that are in the active site. This would allow us to create :ref:`object_definition`.

But wait. Is it really the correct structure? How many chains there are? What is the numeration of residues?

Create *Object deninition*
^^^^^^^^^^^^^^^^^^^^^^^^^^

Lets load another structure. Open file ``first_frame_1cqz.pdb`` downloaded from test data `repository <http://localhost:8001>`_. It is a first frame of the MD simulation and it is en example of how the frame of MD looks like. In order to create :ref:`object_definition` you have to discover folowin things:

#. What is the name of water residue?
#. What are numbers of residues in the active site?
#. What size the active site is?

.. note::

    It is also good idea to open ``.pdb`` file in your favorite text editor and look at residue numbers and names.

Create *Scope deninition*
^^^^^^^^^^^^^^^^^^^^^^^^^^

:ref:`scope_definition` is easy to create. We will use *Convex hull* version so the scope definition could be simply ``backbone``.






