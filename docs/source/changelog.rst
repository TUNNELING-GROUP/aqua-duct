Aqua-Duct changelog
===================

* devel
    * Steady improvement of documentation (including API).
    * Added mario script, WIP.
    * If barber is used for clusterization, apropriate radii are displayed. Buggy.
    * Data dumps routines moved to separare class, will be in separate modeule soon.
    * Tests with netcdf4 for dumping data, WIP.
* 0.3.2
    * Major improvement: new auto_barber based clustering method.
    * Clusterization history displayed as simple ascii tree.
    * AutoBarber min and max cut level options added.
    * Barber moved to separate module.
    * Fixed bug in visualization script; if no molecule is kept do not set style and color.
* 0.3.1
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
* 0.2.26
    * Stage execution time debug messages.
    * Total execution time debug message.
* 0.2.25
    * initial public release
