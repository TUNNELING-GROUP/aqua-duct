Aqua-Duct changelog
===================

* 1.0.2 (29.07.2019)
    * Improvements in GUI configured.
    * Fixed problems with waterfall mode.
* 1.0.0 (21.07.2019) Aqua-Duct version 1
    * Assorted minor fixes and improvements.
* 1.0.0b1 (13.07.2019) beta
    * New Waterfall option for analysis of custom sampled simulations.
    * New Kraken GUI tool for facilitating analysis of results.
    * RAM usage optimization in AutoBarber and other procedures.
    * Center of Object and Center of System calculation.
    * Cluster shapes KDE estimation and visualization.
    * New GUI tool for preparing configuration files.
    * Possibility to join clusters by ID and to sort and renumber cluster IDs.
    * Inlets in seleced clusters can be removed.
    * Order of clustering procedures can be altered by user.
    * AutoBarber can be run for different types of traced molecules separately.
    * New driver Pond allows to calculate pockets and master paths energy profiles.
    * Assorted improvements allowing better handling of passing paths.
    * GREAT speedup of Stage I, II, and III calculations: they run in parallel; IO can be a bottleneck though.
    * Substantial speedup of SinglePaths generation: it runs in parallel and uses fastest routines.
    * Substantial speedup of AutoBarber procedures: it runs in parallel; IO can be a bottleneck though.
    * Improvements in analysis stage. Additional info displayed in tables and added progress bar.
    * portal.py script for calculating sizes of selection(s) using convex hull approximation.
    * Improvements in dir-cache handling.
    * Some speedup of master paths calculations, more to come.
    * Paths, SinglePaths and other objects use less memory.
    * Newest MDAnalysis can be used (ie 0.17, 0.18, 0.19) however it is recommended to stay with 0.16.2.
    * Many other minor improvements and bug fixes.
* 0.5.15 (05.05.2019)
    * Topology residue ID are printed in the analysis files.
* 0.5.14 (26.09.2018)
    * Recommended MDAnalysis is set to >=0.16 and <0.17. Versions >=0.17 are fully supported.
    * Docs update.
    * Various performance improvements and few minor bug fixes.
* 0.5.9 (12.03.2018)
    * Rewritten module for MD data access. Sandwich mode added.
    * Coordinates can be stored in cache directory, in memory or generated on demand.
    * Support for long trajectories.
    * Passing through paths are supported.
    * Improvements in visualization script.
    * Coordinates of residues are calculated as center of geometry.
    * Recommended MDAnalysis is set to >=0.16 and <0.17. Version 0.17 is supported but not recommended.
    * Bug fixes and code cleanup.
* 0.4.0 - 0.4.14 (20.11.2017) unofficial
    * Uses newest MDAnalysis (0.16.2).
    * Steady improvement of documentation (including API).
    * Names of traced molecules are returned in the result file and tables are split appropriately.
    * Tables in the result file are split in regard to Object and Passing paths.
    * Passing through paths are being introduced, WIP.
    * Additional tables in the result file.
    * CRD is enabled as topology/trajectory format.
    * Traced residues are identified by resindices instead of resids; this allows to use weak topologies such as PDB.
    * Removed roman dependency.
    * In addition to histograms approximate (ConvexHull approximation) areas and volumes of the scope and object can be calculated.
    * Bug fixes and reliability fixes.
* 0.3.7 (18.07.2017)
    * Enable XTC trajectory format.
    * Reliability fix in progress bar display.
* 0.3.6 (28.06.2017)
    * AQ can be run for given part of trajectory.
    * Fixed bug in passing options to Barber clustering method.
    * Recursive threshold can be defined as range; no disjoint ranges are supported.
* 0.3.5 (18.04.2017)
    * As for now, the only supported version of MDAnalysis is 0.15.
* 0.3.4 (14.04.2017)
    * Fixed bug in progress bar updating method causing critical error in some specific circumstances.
* 0.3.3 (20.03.2017)
    * AutoBarber default values of maxcut_level and mincut_level changed to True.
    * Improved template configuration file.
    * Number of small improvements in documentation.
* 0.3.2 (24.02.2017)
    * Major improvement: new auto_barber based clustering method.
    * Clustering history displayed as simple ascii tree.
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
