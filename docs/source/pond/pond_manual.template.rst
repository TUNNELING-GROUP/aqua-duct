*Pond* manual
==============

*Pond* application is a driver that uses :mod:`aquaduct` module to perform further analysis of results from *Valve* calculations.

*Pond* can calculate pockets present in the protein and free energy profiles of :ref:`master_paths_manual`.


*Pond* invocation
------------------

Once :mod:`aquaduct` module is installed (see :doc:`../aquaduct_install`) properly on the machine, *Pond* is available as ``pond.py`` command line tool.

Usage
^^^^^

Basic help of *Pond* usage can be displayed by following command::

    pond.py --help

It should display following information::

HELP

Options common with *Valve*
^^^^^^^^^^^^^^^^^^^^^^^^^^^

All options related to Molecular Dynamic simulation data, configuration file, and threads have the same meaning as in *Valve*.

For detailed explanation of the following options see :doc:`../valve/valve_manual`:

* ``-c CONFIG_FILE`` Config file filename. (default: None)
* ``-t THREADS`` Limit Aqua-Duct calculations to given number of threads. (default: None)
* ``--max-frame MAX_FRAME`` Maximal number of frame. (default: None)
* ``--min-frame MIN_FRAME`` Minimal number of frame. (default: None)
* ``--step-frame STEP_FRAME`` Frames step. (default: None)
* ``--sandwich`` Sandwich mode for multiple trajectories. (default: False)
* ``--cache-dir CACHEDIR`` Directory for coordinates caching. (default: None)

*Pond* reference density calculation
------------------------------------

