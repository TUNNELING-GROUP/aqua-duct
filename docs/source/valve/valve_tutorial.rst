*Valve* tutorial
================

This tutorial assumes :mod:`aquaduct` and *Valve* is already installed - see :doc:`../aquaduct_install`. It is also assumed that user is acquainted with :doc:`valve_manual` and *Valve* :doc:`valve_config`.

*Valve* invocation
------------------

Usually *Valve* is run by::

    valve.py

To check if *Valve* is installed and works properly try to issue following commands::

    valve.py --help
    valve.py --version

.. _test_data:

Test data
---------

**Mouse!**

We will use 1 ns MD simulation data of sEH protein (PDBID **1cqz**). This simulation was performed in Amber 14. Necessary files can be found  at `Aqua-Duct home page <http://aquaduct.pl/>`_ in section `download <http://aquaduct.pl/download>`_. Required data is in the `sample data` file.


Inspect your system
-------------------

Before we start any calculations let's have a look at the protein of interest. Start PyMOL and get ``1cqz`` PDB structure (for example by typing in PyMOL command prompt ``fetch 1cqz``).

To setup *Valve* calculations we need to know the active site of the protein. More precisely we need to know IDs of residues that are in the active site. This would allow us to create :ref:`object_definition`.

But wait. Is it really the correct structure? How many chains there are? What is the numeration of residues? How does it compare with the topology file from `sample data`?

Create *Object definition*
^^^^^^^^^^^^^^^^^^^^^^^^^^

Leti's load another structure. Open file ``1cqz_sample_topology.pdb`` (see :ref:`Test data<test_data>`). It is a first frame of the MD simulation and it is an example of how the frame of MD looks like. In order to create :ref:`object_definition` you have to discover following things:

#. What is the name of water molecules?
#. What are numbers of residues in the active site?
#. What size the active site is of?

.. note::

    It is also a good idea to open ``.pdb`` file in your favorite text editor and look at residue numbers and names.

Create *Scope definition*
^^^^^^^^^^^^^^^^^^^^^^^^^^

:ref:`scope_definition` is easy to create. We will use *Convex hull* version so the scope definition could be simply ``backbone``.


Prepare config file
-------------------

*Valve* performs calculations according to the configuration (aka *config*) file.

Lets start from dumping config file template to ``config.txt`` file. Open it in your favorite editor and fill all options.
If you have troubles look at :doc:`valve_config` (and :doc:`valve_manual`).

Things to remember:

#. Provide correct paths to topology and trajectory data.
#. Enter correct :ref:`Object <object_definition>` and :ref:`Scope <scope_definition>` definitions.
#. Make sure visualization is switched on.

Run *Valve*
-----------

Make sure all necessary data is in place. Open terminal, go to your working directory and type in::

    valve.py -c config.txt

Depending on your machine and current load it may take a while (matter of minutes) to complete all calculations.

Visual inspection
^^^^^^^^^^^^^^^^^

In the last stage PyMOL should pop up and *Valve* should start to feed it with visualization data. This would take a moment and if you set up ``save`` option a PyMOL session would be saved. Once it is done *Valve* quits and switches off PyMOL. Now, you can restart it and read saved session.

Clustering
^^^^^^^^^^^^^^

Improve clustering of Inlets. See :doc:`valve_config` for more hints on available clustering options.

Analysis tables
^^^^^^^^^^^^^^^

Open ``5_analysis_results.txt`` file and look at summaries and tables. See also :doc:`valve_manual`.

Feedback
--------

Give us your opinion. Send your questions, inquires, anything to developer(s): `info@aquaduct.pl <info@aquaduct.pl>`_.
There are couple of questions that might be useful to form your opinion.

#. What do you like in *Valve* and *Aqua-Duct*?
#. What do you do not like in *Valve* or *Aqua-Duct*?
#. What is missing?
#. Do you find it useful?
