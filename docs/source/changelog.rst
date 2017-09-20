Aqua-Duct changelog
===================

* devel
    * Uses newest MDAnalysis (0.16.2).
    * Steady improvement of documentation (including API).
    * Added mario script, WIP.
    * If barber is used for clusterization, appropriate radii are displayed. Buggy.
    * Data dumps routines moved to separate module.
    * Tests with netcdf4 for dumping data, WIP.
    * SmartRanges moved to helpers module.
    * Names of traced molecules are returned in the result file and tables are split appropriately.
    * Tables in the result file are split in regard to Object and Passing paths.
    * Small bug in reporting progress in AutoBarber preparation fixed.
    * Passing through paths are being introduced, WIP.
    * Small bug in linearization functions fixed.
    * CRD is enabled as topology/trajectory format.
    * Traced residues are identified by resindices instead of resids; this allows to use weak topologies such as PDB.
    * Stage II performance improved if clear_in_object_info is not used (default).
* 0.3.7 (18.07.2017)
    * Enable XTC trajectory format.
    * Reliability fix in progress bar display.
* 0.3.6 (28.06.2017)
    * AQ can be run for given part of trajectory.
    * Fixed bug in passing options to Barber clusterization method.
    * Recursive threshold can be defined as range; no disjoint ranges are supported.
* 0.3.5 (18.04.2017)
    * As for now, the only supported version of MDAnalysis is 0.15.
* 0.3.4 (14.04.2017)
    * Fixed bug in progress bar updating method causing critical error in some specific circumstances.
* 0.3.3 (20.03.2017)
    * AutoBarber default values of maxcut_level and mincut_level changed to True.
    * Improved template configuration file.
    * Number of small improvements in documentaion.
* 0.3.2 (24.02.2017)
    * Major improvement: new auto_barber based clustering method.
    * Clusterization history displayed as simple ascii tree.
    * AutoBarber min and max cut level options added.
    * Barber moved to separate module.
    * Fixed bug in visualization script; if no molecule is kept do not set style and color.
* 0.3.1 (04.02.2017)
    * AutoBarber tovdw option.
    * AutoBarber minimal and maximal cut options.
    * Fixed bug in AutoBarber: some areas were sometimes not cut.
    * Documentation improvements.
    * Valve driver simplified. Most of the functionality moved to separate module.
    * Option for single precision storage.
    * Added Savitzky-Golay smoothing; AQ requires SciPy >= 0.14 now.
    * Improved sorting of CTypes.
    * Raw and Separate paths uses SmartRanges. This allowed for excellent performance improvement of Separate paths calculation.
    * Default display of molecule changed to silver cartoon.
    * Object shape displayed in orange.
    * Fixed several small bugs.
* 0.2.26 (21.01.2017)
    * Stage execution time debug messages.
    * Total execution time debug message.
* 0.2.25 (18.01.2017)
    * initial public release
