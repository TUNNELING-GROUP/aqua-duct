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

    usage: pond.py [-h] [-c CONFIG_FILE] [-t THREADS] [-r RESULTS_DIR]
                   [--max-frame MAX_FRAME] [--min-frame MIN_FRAME]
                   [--step-frame STEP_FRAME] [--raw] [--raw-master]
                   [--raw-discard-singletons RAW_SINGL] [--sandwich]
                   [--cache-dir CACHEDIR] [--window-full] [--windows WINDOWS]
                   [--wsize WSIZE] [--reference REF]
                   [--reference-radius REF_RADIUS] [--reference-mol REF_MOL]
                   [--temperature TEMP] [--gsize GRID_SIZE] [--pockets]
                   [--hotspots] [--master-radius MASTER_RADIUS]
                   [--master-ctypes MASTER_CTYPES]
    
    What have I got in my pocket?
    
    optional arguments:
      -h, --help            show this help message and exit
      -c CONFIG_FILE        Config file filename. (default: None)
      -t THREADS            Limit Aqua-Duct calculations to given number of
                            threads. (default: None)
      -r RESULTS_DIR        Path to results directory (default: )
      --max-frame MAX_FRAME
                            Maximal number of frame. (default: None)
      --min-frame MIN_FRAME
                            Minimal number of frame. (default: None)
      --step-frame STEP_FRAME
                            Frames step. (default: None)
      --raw                 Use raw data from paths instead of single paths.
                            (default: False)
      --raw-master          Use raw data from paths instead of single paths, only
                            in master paths calculations. (default: False)
      --raw-discard-singletons RAW_SINGL
                            Discard short scope only segments from raw data.
                            (default: 1)
      --sandwich            Sandwich mode for multiple trajectories. (default:
                            False)
      --cache-dir CACHEDIR  Directory for coordinates caching. (default: None)
      --window-full         Return full window if windows is used. (default:
                            False)
      --windows WINDOWS     Number of windows to calculate. (default: 1)
      --wsize WSIZE         Size of window in frames. (default: None)
      --reference REF       Selection of reference in the first frame of
                            trajectory. (default: None)
      --reference-radius REF_RADIUS
                            Radius of reference. (default: 2.0)
      --reference-mol REF_MOL
                            Selection of reference molecules. (default: resname
                            WAT)
      --temperature TEMP    Simulation temperature. (default: 300.0)
      --gsize GRID_SIZE     Size of grid's cells. (default: 1.0)
      --pockets             Calculate pockets. (default: False)
      --hotspots            Calculates hotspots if pockets are calculated.
                            (default: False)
      --master-radius MASTER_RADIUS
                            Calculate profiles for master paths with given radius.
                            (default: None)
      --master-ctypes MASTER_CTYPES
                            Limit calculations to given ctypes. (default: None)

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

