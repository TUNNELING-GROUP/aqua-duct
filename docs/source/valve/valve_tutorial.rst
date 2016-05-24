*Valve* tutorial
================

This is a tentative *Valve* manual. Created for the sake of Aqueduct training we have today. Eventually, it will be rewritten to the official version.

This tutorial assumes :mod:`aqueduct` and *Valve* is already installed - see :doc:`../aqueduct_install`. It is also assumed that user is acquainted with :doc:`valve_manual` and *Valve* :doc:`valve_config`.


*Valve* invocation
------------------

Usually *Valve* is run by::

    valve.py

Due to specific setup in our laboratory *Valve* has to be run through simple wrapper script::

    valve_run

Additionally, to speed up all calculations it is assumed that *Valve* is run with ``--max-frame 1000`` option::

    valve_run --max-frame 1000

To check is *Valve* is installed and works properly try to issue following commands::

    valve_run --help
    valve_run --version

Test data
---------

**Mouse!**

We will use 10ns Amber MD simulation data of sEH protein (PDBID **1cqz**). Necessary files can be downloaded `here <http://localhost:8001>`_:

* Go to download server.
* Go to ``1cqz`` directory.
* Download all files and save them in sane location on your machine. Please note, that ``.nc`` file is ca. 3.5 GB so it may take a while to download it.

Inspect your system
-------------------

Before we start any calculations lets have a look at the protein of interest. Start *PyMOL* and get ``1cqz`` PDB structure (for example by typing in *PyMOL* command prompt ``fetch 1cqa``).

To setup *Valve* calculations we need to know active site of the protein. More precisely we need to know IDs or residues that are in the active site. This would allow us to create :ref:`object_definition`.

But wait. Is it really the correct structure? How many chains there are? What is the numeration of residues?

Create *Object definition*
^^^^^^^^^^^^^^^^^^^^^^^^^^

Lets load another structure. Open file ``first_frame_1cqz.pdb`` downloaded from test data `repository <http://localhost:8001>`_. It is a first frame of the MD simulation and it is en example of how the frame of MD looks like. In order to create :ref:`object_definition` you have to discover following things:

#. What is the name of water residue?
#. What are numbers of residues in the active site?
#. What size the active site is?

.. note::

    It is also good idea to open ``.pdb`` file in your favorite text editor and look at residue numbers and names.

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
#. Provide file name of result in analysis section, for example ``results.txt`` (for future reference).
#. Make sure visualization is switched on and ``save`` option points to session file name (``.pse``)

Run *Valve*
-----------

Make sure all necessary data is in place. Open terminal, go to your working directory and type in::

    valve_run --max-frame 1000 -c config.txt

Depending on your machine and current load it may take a while (matter of minutes) to complete all calculations.

Visual inspection
^^^^^^^^^^^^^^^^^

In the last stage *PyMOL* should pop up and *Valve* should start to feed it with visualization data. This would take a moment and if you set up ``save`` option a *PyMOL* session would be saved. Once it is done *Valve* quits and switches off *PyMOL*. Now, you can restart it and read saved session.

Analysis tables
^^^^^^^^^^^^^^^

Open ``results.txt`` file and look at summaries and tables. See also :doc:`valve_manual`.

Feedback
--------

Give us your opinion. Send your questions, inquires, anything to developer(s): `Tomasz Magdziarz <t.magdziarz@tunnelinggroup.pl>`_.
This are couple of questions that might be useful to form your opinion.

#. What do you like in *Valve* and *Aqueduct*?
#. What do you do not like in *Valve* or *Aqueduct*?
#. What is missing?
#. Do you find it useful?

