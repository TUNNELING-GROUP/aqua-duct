Aqua-Duct changelog
===================

* devel
    * Added mario script, WIP.
    * Fixed bug in visualization script; if no molecule is kept do not set style and color.
    * Barber moved to separate module.
    * Bug identified, fix pending: reclusterization is not run in a recursive manner.
    * Major improvement: new auto_barber based clustering method.
    * If barber is used for clusterization, apropriate radii are displayed.
    * Clusterization history displayed as simple ascii tree.
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
