Configuration file options
==========================

Valve configuration file is a simple and plain text file. It has similar structure as INI files commonly used in one of the popular operating systems and is compliant with Python module :mod:`ConfigParser`.

Configuration file comprises of several *sections*. They can be grouped into three categories. Names of sections are in **bold** text.

#. Global settings:
    * **global**
#. Stages options:
    #. **traceable_residues**
    #. **raw_paths**
    #. **separate_paths**
    #. **inlets_clustering**
    #. **analysis**
    #. **visualize**
#. Methods options:
    * **smooth**
    * **clustering**
    * **reclusteriation**

Section **global**
------------------

This section allows settings of trajectory data and is reserved for other future global options.

=========== =============   ==========================================================================
Option      Default value   Description
=========== =============   ==========================================================================
top         None            Path to topology file. Aqua-Duct supports PDB, PRMTOP, PFS topology files.
trj         None            Path to trajectory file. Aqua-Duct supports NC and DCD trajectory files.
twoway      True            Try to use *two-way* scanning in the stage II.
sandwich    False           If set ``True`` trajectories are read as layers.
max_frame   None            Maximal number of frame to be read from trajectory data. If set ``None``
                            trajectory is read to the last possible frame.
min_frame   0               Minimal number of frame to be read from trajectory data.
step_frame  1               Step used in reading trajectory. Default value of 1 stands for reading
                            every frame. If it is greater than 1, only every step-value frame is
                            read.
sps         True            Try to store data in single precission storage.
cache_dir   None            Allows to set path to the directory for cache data.
cache_mem   False           If set ``True``, all data will be cached in RAM.
=========== =============   ==========================================================================

Option **trj** can be used to provide list of trajectory files separated by standard path separator '``:``' on POSIX platforms and '``;``' on Windows - see :obj:`os.pathsep`.

.. note::

    Options **top** and **trj** are mandatory.


.. note::

    Options **min_frame**, **max_frame**, and **step_frame** can be used to limit calculations to a specific part of trajectory. For example, in order to to run calculations for 1000 frames starting from frame 5000 use the following options::

        min_frame = 4999
        max_frame = 5999


    To run calculations for every 5th frame use::

        step_frame = 5


.. _sandwich_option:

Sandwich
^^^^^^^^

Trajectory data can be provided as several files. By default these files are processed in sequential manner making one long trajectory. If option **sandwich** is used trajectory files are read as layers. For each layer, search for traceable residues is done separately (stage I and II) but processing and analysis (stage III, IV, V, and VI) are done for all paths simultaneously. Usage of **sandwich** option is further referenced as *sandwich* mode.


Cache
^^^^^

Storage of coordinates for all paths for very long MD trajectories requires huge amount of RAM. User can decide whether :mod:`aquaduct` should store coordinates in memory or in separated directory. Option **cache-mem** instruct *Valve* to store coordinates in RAM; **cache-dir** stores coordinates in selected directory. If neither of both options is selected, coordinates are calculated on demand.

.. note::

    If no cache is used (memory or dir) :ref:`master_paths_manual` cannot be calculated.


Single precision storage
^^^^^^^^^^^^^^^^^^^^^^^^

Most of the calculation in *Valve* is performed by NumPy. By default, NumPy uses double precision floats.
*Valve* does not change this behavior but has special option **sps** which forces to store all data (both internal data stored in RAM and on the disk) in single precision. This spares a lot of RAM and is recommended when you perform calculation for long trajectories and your amount of RAM is limited.



Common settings of stage sections
---------------------------------

Stages 1-4 which perform calculations have some common options allowing for execution control and saving/loading data.

.. tabularcolumns:: |p{1.0cm}|p{2.5cm}|p{11.1cm}|

========    =================   ===================================================================
Option      Default value       Description
========    =================   ===================================================================
execute     runonce             Option controls stage execution. It can have one of three possible
                                values: ``run``, ``runonce``, and ``skip``. If it is set to ``run``
                                calculations are always performed and if **dump** is set, dump file
                                is saved. If it is set to ``runonce``, calculations are performed
                                when there is no dump file specified by **dump** option. If it is
                                present, calculations are skipped and data is loaded from the file.
                                If it is set to ``skip``, calculations are skipped and if **dump**
                                is set, data is loaded from the file.
dump        [dump file name]    File name of dump data. It is used to save results of calculations
                                or to load previously calculated data - this depends on **execute**
                                option. Default value of this option depends on the stage and for
                                stages 1 to 4 is one of the following (listed in order):

                                * 1_traceable_residues_data.dump
                                * 2_raw_paths_data.dump
                                * 3_separate_paths_data.dump
                                * 4_inlets_clustering_data.dump
========    =================   ===================================================================

Stages 5-6 also use **execute** option, however, since they do not perform calculations `per se`, instead of **dump** option, they use **save**.

.. tabularcolumns:: |p{1.0cm}|p{2.5cm}|p{11.1cm}|

========    =================   ===================================================================
Option      Default value       Description
========    =================   ===================================================================
execute     run                 Option controls stage execution. It can have one of three possible
                                values: ``run``, ``runonce``, and ``skip``. If it is set to ``run``
                                or ``runonce``, stage is executed and results are saved according to
                                the **save** option. If it is set to ``skip``, stage is skipped.
save        [save file name]    File name for saving results. Default value of this option depends
                                on the stage and for stages 1 to 4 is one of the following
                                (listed in order):

                                * 5_analysis_results.txt & 5_analysis_results.txt.csv
                                * 6_visualize_results.py & 6_visualize_results.tar.gz

                                Stage 5 saves ``.txt`` file with analysis of results and, if
                                requested, it saves additional ``.csv`` with various counts of
                                traced molecules.

                                Stage 6 can save results in two different ways:

                                #. As Python script - extension ``.py`` plus companion archive
                                   ``.tar.gz``,
                                #. As PyMOL session - extension ``.pse``.
========    =================   ===================================================================


Stage **traceable_residues**
----------------------------

.. _inflate_options:


.. tabularcolumns:: |p{3.6cm}|p{2.5cm}|p{8.5cm}|

=========================   ==============  ================================================================
Option                      Default value   Description
=========================   ==============  ================================================================
scope                       None            Definition of *Scope* of interest. See also
                                            :ref:`scope_definition`.
scope_convexhull            True            Flag to set if *Scope* is direct or convex hull definition.
scope_everyframe            False           Flag to set *Scope* evaluation mode. If set ``True`` *Scope* is
                                            evaluated in every frame. This makes sense if the definition is
                                            complex and depends on distances between molecular entities.
scope_convexhull_inflate    None            Increase (or if negative - decrease) size of the scope convex
                                            hull.
object                      None            Definition of *Object* of interest. See also
                                            :ref:`object_definition`.
add_passing                 None            Definition of molecules that should be added to traced molecules
                                            even if they were not present in *Object*.
=========================   ==============  ================================================================


.. note::

    Options **scope** and **object** are mandatory.


Stage **raw_paths**
-------------------

This stage also requires definition of the *Scope* and *Object*. If appropriate settings are not given, settings from the previous stage are used.

.. tabularcolumns:: |p{3.6cm}|p{2.5cm}|p{8.5cm}|

=========================   ==============  ================================================================
Option                      Default value   Description
=========================   ==============  ================================================================
scope                       None            Definition of *Scope* of interest. See also
                                            :ref:`scope_definition`. If ``None``, value from previous stage
                                            is used.
scope_convexhull            None            Flag to set if the *Scope* is direct or convex hull definition.
scope_everyframe            False           Flag to set *Scope* evaluation mode. If set ``True`` *Scope* is
                                            evaluated in every frame. This makes sense if the definition is
                                            complex and depends on distances between molecular entities.
                                            If ``None``, value from previous stage is used.
scope_convexhull_inflate    None            Increase (or if negative - decrease) size of the scope convex
                                            hull. If ``None``, value from previous stage is used.
object                      None            Definition of *Object* of interest. See also
                                            :ref:`object_definition`. If ``None``, value from previous
                                            stage is used
clear_in_object_info        False           If it is set to ``True``, information on occupation of *Object*
                                            site by traceable residues calculated in the previous stage is
                                            cleared and has to be recalculated. This is useful if
                                            definition of *Object* was changed.
discard_singletons          1               If ``> 0``, discards paths of given lenght.
discard_empty_paths         True            If set to ``True``, empty paths are discarded.
=========================   ==============  ================================================================

.. _separate_paths_options:

Stage **separate_paths**
------------------------

.. tabularcolumns:: |p{4.0cm}|p{2.5cm}|p{8.1cm}|

========================    ==============  ================================================================
Option                      Default value   Description
========================    ==============  ================================================================
discard_empty_paths         True            If set to ``True``, empty paths are discarded.
sort_by_id                  True            If set to ``True``, separate paths are sorted by ID. Otherwise
                                            they are sorted in order of appearance.
discard_short_paths         20              This option allows to discard paths which are shorter than the
                                            threshold, which is defined as the total number of frames.
discard_short_object        2.0             This option allows to discard paths whose objects are shorter
                                            than the threshold, which is defined as total length in metric
                                            units.
discard_short_logic         or              If both ``discard_short_paths`` and ``discard_short_object``
                                            options are used, this option allows to set combination logic.
                                            If it is set to ``or``, a path is discarded if any of discard
                                            criterion is met. If it is set ``and``, both criteria have to
                                            be met to discard path.
auto_barber                 None            This option allows to select molecular entity used in Auto
                                            Barber procedure. See also :ref:`auto_barber_procedure` and
                                            :meth:`~aquaduct.traj.paths.GenericPaths.barber_with_spheres`.
auto_barber_mincut          None            Minimal radius of spheres used in Auto Barber. If a sphere has
                                            radius smaller than this value, it is not used in AutoBarber
                                            procedure. This option can be switched off by setting it to
                                            ``None``.
auto_barber_maxcut          2.8             Maximal radius of spheres used in Auto Barber. If a sphere has
                                            radius greater than this value, it is not used in AutoBarber
                                            procedure. This option can be switched off by setting it to
                                            ``None``.
auto_barber_mincut_level    True            If set ``True``, spheres of radius smaller than **mincut** are
                                            resized to **mincut** value.
auto_barber_maxcut_level    True            If set ``True``, spheres of radius greater than **maxcut** are
                                            resized to **maxcut** value.
auto_barber_tovdw           True            If set ``True``, cutting of spheres is corrected by decreasing
                                            its radius by VdW radius of the closest atom.
allow_passing_paths         False           If set ``True``, paths that do not enter the object are detected
                                            and added to the rest of paths as 'passing' paths.
separate_barber             True            Apply AutoBarber for each type of traced molecules separately.
calculate_coo               False           If set ``True``, center of the object is calculated.
========================    ==============  ================================================================


Stage **inlets_clustering**
-------------------------------

.. tabularcolumns:: |p{5.0cm}|p{2.5cm}|p{7.1cm}|

==================================  ==============  ================================================================
Option                              Default value   Description
==================================  ==============  ================================================================
recluster_outliers                  False           If set to ``True``, reclustering of outliers is executed
                                                    according to the method defined in **reclustering** section.
detect_outliers                     False           If set, detection of outliers is executed. It could be set as a
                                                    floating point distance threshold or set to ``Auto``. See
                                                    :ref:`clustering_of_inlets` for more details.
singletons_outliers                 False           Maximal size of cluster to be considered as outliers. If set to
                                                    number > 0 clusters of that size are removed and their objects
                                                    are moved to outliers. See :ref:`clustering_of_inlets` for
                                                    more details.
max_level                           5               Maximal number of recursive clustering levels.
create_master_paths                 False           If set to ``True``, master paths are created (fast CPU and big
                                                    RAM recommended; 50k frames long simulation may need ca 20GB of
                                                    memory)
master_paths_amount                 None            Allows to limit number of single paths used for master paths
                                                    calculations.
                                                    If it is a number
                                                    in range ``(0,1)``, then it is interpreted as a percent number
                                                    of paths to be used. It is is a integer number ``>= 1``
                                                    it is an absoulte number of used paths.
separate_master                     False           If set to ``True``, master mapths are generated for each of
                                                    traced molecules' types separately.
separate_master_all                 True            If **separate_master** is used and this option is set ``True``,
                                                    master paths will be also generated for all traced molecules`
                                                    types as it is done when **separate_master** is set ``False``.
exclude_passing_in_clustering       True            If set to ``True``, passing paths are not clustered with normal
                                                    paths.
add_passing_to_clusters             None            Allows to run procedure for adding passing paths inlets to
                                                    clusters with Auto Barber method. To enable this the option
                                                    should be set to molecular entity that will be used by Auto
                                                    Barber.
renumber_clusters                   False           If set ``True``, clusters have consecutive numbers starting from
                                                    1 (or 0 if outliers are present) starting from the bigest
                                                    cluster.
join_clusters                       None            This option allows to join selected clusters. Clusters' IDs
                                                    joined with ``+`` character lists clusters to be joined
                                                    together. Several such blocks separated by space can be used.
                                                    For example, if set to ``1+3+4 5+6`` clusters 1, 3, and 4 will
                                                    be joined in one cluster and cluster 5, and 6 will be also
                                                    joined in another cluster.
inlets_center                       cos             Allows to choose center of inlets points. This central point is
                                                    further used as a reference point in calculations of clusters'
                                                    areas and contours. If set `cos` center of system calculated as
                                                    average center of the scope area is used. Alternatively it can
                                                    be set to `coo`, then center of the object area is used.
clustering_order                    old-school      Allow to change order of clustering steps.

                                                    * **old-school** order - clustering, reclustering, removing singletons, adding passing paths to clusters, join clusters, renumbering clusters, removing inlets in specified_clusters
                                                    * **aquarius** order - clustering, joining clusters, renumbering clusters, outliers detection, reclustering, removing singletons, adding passing paths to clusters, removing inlets in specified_clusters
==================================  ==============  ================================================================

Stage **analysis**
------------------

.. tabularcolumns:: |p{4.5cm}|p{2.5cm}|p{7.6cm}|

==============================  ==============  ================================================================
Option                          Default value   Description
==============================  ==============  ================================================================
dump_config                     True            If set to ``True``, configuration options, as seen by Valve, are
                                                added to the head of results.
calculate_scope_object_size     False           If set to ``True``, volumes and areas of object and scope
                                                approximated by convex hulls are calculated for each of the
                                                analyzed frames and saved in output CSV file.
scope_chull                     None            Scope convex hull definition used in calculating volume and
                                                area.
scope_chull_inflate             None            Increase (or if negative - decrease) size of the scope convex
                                                hull.
object_chull                    None            Object convex hull definition used in calculating volume and
                                                area.
cluster_area                    True            If set ``True``, clusters' areas are estimated with kernel
                                                density estimation method (KDE).
cluster_area_precision          20              Precision of KDE method in clusters' areas estimation method.
                                                This options controls number of grid points per one
                                                square A as used in KDE. Higher values means better precision.
                                                Number of points can be calculated as $P^{2/3}$.
cluster_area_expand             2               Space occupied by clusters' points can be expanded before KDE
                                                calculation. This option controls amount of A by which the
                                                cluster space is expanded.
                                                Average amount of expansion can be calcualted as $E^{2/3}$.
==============================  ==============  ================================================================


Stage **visualize**
-------------------

.. tabularcolumns:: |p{4.0cm}|p{2.5cm}|p{8.1cm}|

.. table::
    :class: longtable

    ==========================  ================    ==========================================================================================
    Option                      Default value       Description
    ==========================  ================    ==========================================================================================
    split_by_type               False               If ``True``, visulaizations of all objects that correspond to particular types of
                                                    traced molecules are splited. Appropriate molecule name is added to created objects.
    retain_all_types            False               If ``True`` and *split_by_type* option is ``True`` as well, standard visualization for
                                                    all types of traced moelcules will be generated as well.

    all_paths_raw               False               If ``True``, produces one object in PyMOL that holds all paths
                                                    visualized by raw coordinates.
    all_paths_smooth            False               If ``True``, produces one object in PyMOL that holds all paths
                                                    visualized by smooth coordinates.
    all_paths_split             False               If is set ``True``, objects produced by **all_paths_raw** and
                                                    **all_paths_smooth** are split into Incoming, Object, and
                                                    Outgoing parts and visualized as three different objects.
    all_paths_raw_io            False               If set ``True``, arrows pointing beginning and end of paths are
                                                    displayed oriented accordingly to raw paths orientation.
    all_paths_smooth_io         False               If set ``True``, arrows pointing beginning and end of paths are
                                                    displayed oriented accordingly to smooth paths orientation.
    all_paths_amount            None                Allows to limit number of visualised paths. If it is a number
                                                    in range ``(0,1)``, then it is interpreted as a percent number
                                                    of paths to be visualized. It is is a integer number ``>= 1``
                                                    it is total number of all_paths visualized.
    simply_smooths              RecursiveVector     Option indicates linear simplification method to be used in
                                                    plotting smooth paths. Simplification removes points which do
                                                    not (or almost do not) change the shape of smooth path.
                                                    Possible choices are:

                                                    * ``RecursiveVector`` (:class:`~aquaduct.geom.traces.LinearizeRecursiveVector`),
                                                    * ``HobbitVector`` (:class:`~aquaduct.geom.traces.LinearizeHobbitVector`),
                                                    * ``OneWayVector`` (:class:`~aquaduct.geom.traces.LinearizeOneWayVector`),
                                                    * ``RecursiveTriangle`` (:class:`~aquaduct.geom.traces.LinearizeRecursiveTriangle`),
                                                    * ``HobbitTriangle`` (:class:`~aquaduct.geom.traces.LinearizeHobbitTriangle`),
                                                    * ``OneWayTriangle`` (:class:`~aquaduct.geom.traces.LinearizeOneWayTriangle`).

                                                    Optionally name of the method can be followed by a threshold
                                                    value in parentheses, i.e. ``RecursiveVector(0.05)``. For sane
                                                    values of thresholds see appropriate documentation of each method.
                                                    Default values work well. This option is not case sensitive.
                                                    It is recommended to use default method or ``HobbitVector`` method.
    paths_raw                   False               If set ``True``, raw paths are displayed as separate objects or as
                                                    one object with states corresponding to number of path.
    paths_smooth                False               If set ``True``, smooth paths are displayed as separate objects or
                                                    as one object with states corresponding to number of path.
    paths_raw_io                False               If set ``True``, arrows indicating beginning and end of paths,
                                                    oriented accordingly to raw paths, are displayed as separate
                                                    objects or as one object with states corresponding to number
                                                    of paths.
    paths_smooth_io             False               If set ``True``, arrows indicating beginning and end of paths,
                                                    oriented accordingly to smooth paths, are displayed as separate
                                                    objects or as one object with states corresponding to number
                                                    of paths.
    paths_states                False               If set ``True``, objects displayed by **paths_raw**,
                                                    **paths_smooth**, **paths_raw_io**, and **paths_smooth_io** are
                                                    displayed as one object with states corresponding to number of
                                                    paths. Otherwise they are displayed as separate objects.
    ctypes_raw                  False               Displays raw paths in a similar manner as non split
                                                    **all_paths_raw** but each cluster type is displayed as
                                                    a separate object.
    ctypes_smooth               False               Displays smooth paths in a similar manner as non split
                                                    **all_paths_smooth** but each cluster type is displayed as
                                                    a separate object.
    ctypes_amount               None                Allows to limit number of visualised ctypes. If it is a number
                                                    in range ``(0,1)``, then it is interpreted as percent number
                                                    of ctypes to be visualized. It is is a integer number ``>= 1``,
                                                    it is total number of visualized ctypes.
    inlets_clusters             False               If set ``True``, clusters of inlets are visualized.
    inlets_clusters_amount      None                Allows to limit number of visualised inlets. If it is a number
                                                    in range ``(0,1)`` then it is interpreted as percent number
                                                    of inlets to be visualized. It is is a integer number ``>= 1``
                                                    it is total number of visualized inlets.
    show_molecule               False               If set to selection of some molecular object in the system,
                                                    for example to ``protein``, this object is displayed.
    show_molecule_frames        0                   Allows to indicate which frames of object defined by
                                                    **show_molecule** should be displayed. It is possible to set
                                                    several frames. In that case frames would be displayed as
                                                    states.
    show_scope_chull            False               If set to selection of some molecular object in the system,
                                                    for example to ``protein``, convex hull of this object is
                                                    displayed.
    show_scope_chull_inflate    None                Increase (or if negative decrease) size of the scope convex
                                                    hull.
    show_scope_chull_frames     0                   Allows to indicate for which frames of object defined by
                                                    **show_chull** convex hull should be displayed. It is possible
                                                    to set several frames. In that case frames would be displayed
                                                    as states.
    show_object_chull           False               If set to selection of some molecular object in the system,
                                                    convex hull of this object is displayed. This works exacly the
                                                    same way as **show_chull** but is meant to mark object shape.
                                                    It can be achieved by using `name * and` molecular object
                                                    definition plus some spatial constrains, for example those
                                                    used in object definition.
    show_object_chull_frames    0                   Allows to indicate for which frames of object defined by
                                                    **show_object** convex hull should be displayed. It is possible
                                                    to set several frames. In that case frames would be displayed
                                                    as states.
    cluster_area                True                If set ``True``, clusters' areas are estimated with kernel
                                                    density estimation method (KDE) and plotted as countour.
    cluster_area_precision      20                  Precision of KDE method in clusters' areas estimation method.
                                                    This options controls number of grid points per one
                                                    square A as used in KDE. Higher values means better precision.
                                                    Number of points can be calculated as $P^{2/3}$.
    cluster_area_expand         2                   Space occupied by clusters' points can be expanded before KDE
                                                    calculation. This option controls amount of A by which the
                                                    cluster space is expanded.
                                                    Average amount of expansion can be calcualted as $E^{2/3}$.
    ==========================  ================    ==========================================================================================


.. note::

    Possibly due to limitations of :mod:`MDAnalysis` only whole molecules can be displayed. If **show_molecule** is set to ``backbone`` complete protein will be displayed anyway. This may change in future version of :mod:`MDAnalysis` and or :mod:`aquaduct`.

.. note::

    If several frames are selected, they are displayed as states which may interfere with other PyMOL objects displayed with several states.

.. note::

    If several states are displayed, protein tertiary structure data might be lost. This seems to be limitation of either :mod:`MDAnalysis` or PyMOL.

.. _clustering_options:

Clustering sections
-----------------------

Default section for definition of clustering method is named **clustering** and default section for reclustering method definition is named **reclustering**. All clustering sections shares some common options. Other options depends on the method.

.. tabularcolumns:: |p{3.5cm}|p{2.5cm}|p{8.6cm}|

=========================   =============== ================================================================
Option                      Default value   Description
=========================   =============== ================================================================
method                      barber or       Name of clustering method. It has to be one of the
                            dbscan          following: barber, dbscan, affprop, meanshift, birch, kmeans.
                                            Default value depends whether it is **clustering** section
                                            (barber) or **reclustering** section (dbscan).
recursive_clustering        clustering      If set to name of some section that holds clustering
                            or None         method settings, this method is called in the next
                                            recursion of clusteriation. Default value for
                                            **reclustering** is None.
recursive_threshold         None            Allows to set threshold that excludes clusters of certain
                                            size from reclustering. Value of this option comprises of
                                            `operator` and `value`. Operator can be one of the following:
                                            >, >=, <=, <. Value has to be expressed as floating number and
                                            it have to be in the range of 0 to 1. One can use several
                                            definitions separated by a space character.
                                            Only clusters of size complying with all thresholds definitions
                                            are submitted to reclustering.
=========================   =============== ================================================================

.. _clustering_methods:

barber
^^^^^^

.. _clustering_barber:


Clustering by **barber** method bases on :ref:`auto_barber_procedure` procedure. For each inlets a sphere is constructed according to Auto Barber **separate_paths** stage settings or according to parameters given in clustering section. Next, inlets that form coherent clouds of mutually intersecting spheres are grouped into clusters. Method **barber** supports the same settings as Auto Barber settings:


.. tabularcolumns:: |p{4.0cm}|p{2.5cm}|p{8.1cm}|

========================    ==============  ================================================================
Option                      Value type      Description
========================    ==============  ================================================================
auto_barber                 str             This option allows to select molecular entity used in Auto
                                            Barber procedure. See also :ref:`auto_barber_procedure` and
                                            :meth:`~aquaduct.traj.paths.GenericPaths.barber_with_spheres`.
auto_barber_mincut          float           Minimal radius of spheres used in Auto Barber. If a sphere has
                                            radius smaller than this value, it is not used to cut. This
                                            option can be switched off by setting it to ``None``.
auto_barber_maxcut          float           Maximal radius of spheres used in Auto Barber. If a sphere has
                                            radius greater than this value, it is not used to cut. This
                                            option can be switched off by setting it to ``None``.
auto_barber_mincut_level    bool            If set ``True``, spheres of radius less than **mincut** are
                                            resized to **mincut** value.
auto_barber_maxcut_level    bool            If set ``True``, spheres of radius greater than **maxcut** are
                                            resized to **maxcut** value.
auto_barber_tovdw           bool            If set ``True``, cutting of spheres is corrected by decreasing
                                            its radius by VdW radius of the closest atom.
========================    ==============  ================================================================

dbscan
^^^^^^

For detailed description look at :class:`sklearn.cluster.DBSCAN` documentation. The following table summarizes options available in `Valve` and is a copy of original documentation.

.. tabularcolumns:: |p{2.5cm}|p{2.5cm}|p{9.6cm}|

=========================   =============== ================================================================
Option                      Value type      Description
=========================   =============== ================================================================
eps                         float           The maximum distance between two samples for them to be
                                            considered as in the same neighborhood.
min_samples                 int             The number of samples (or total weight) in a neighborhood for
                                            a point to be considered as a core point. This includes the
                                            point itself.
metric                      str             The metric to use when calculating distance between instances
                                            in a feature array. Can be one of the following:

                                            * ``euclidean``,
                                            * ``cityblock``,
                                            * ``cosine``,
                                            * ``manhattan``.
algorithm                   str             The algorithm to be used by the NearestNeighbors module to
                                            compute pointwise distances and find nearest neighbors.
                                            Can be one of the following:

                                            * ``auto``,
                                            * ``ball_tree``,
                                            * ``kd_tree``,
                                            * ``brute``.
leaf_size                   int             Leaf size passed to BallTree or cKDTree.
=========================   =============== ================================================================

affprop
^^^^^^^

For detailed description look at :class:`~sklearn.cluster.AffinityPropagation` documentation. The following table summarizes options available in `Valve` and is a copy of original documentation.

.. tabularcolumns:: |p{2.5cm}|p{2.5cm}|p{9.6cm}|

=========================   =============== ================================================================
Option                      Value type      Description
=========================   =============== ================================================================
damping                     float           Damping factor between 0.5 and 1.
convergence_iter            int             Number of iterations with no change in the number of estimated
                                            clusters that stops the convergence.
max_iter                    int             Maximum number of iterations.
preference                  float           Points with larger values of preferences are more likely to be
                                            chosen as exemplars.
=========================   =============== ================================================================

meanshift
^^^^^^^^^

For detailed description look at :class:`~sklearn.cluster.MeanShift` documentation. Following table summarized options available in `Valve` and is a copy of original documentation.

.. tabularcolumns:: |p{2.5cm}|p{2.5cm}|p{9.6cm}|

=========================   =============== ================================================================
Option                      Value type      Description
=========================   =============== ================================================================
bandwidth                   Auto or float   Bandwidth used in the RBF kernel. If ``Auto`` or ``None``
                                            automatic method for bandwidth estimation is used. See
                                            :func:`~sklearn.cluster.estimate_bandwidth`.
cluster_all                 bool            If true, then all points are clustered, even those orphans that
                                            are not within any kernel.
bin_seeding                 bool            If true, initial kernel locations are not locations of all
                                            points, but rather the location of the discretized version of
                                            points, where points are binned onto a grid whose coarseness
                                            corresponds to the bandwidth.
min_bin_freq                int             To speed up the algorithm, accept only those bins with at least
                                            min_bin_freq points as seeds. If not defined, set to 1.
=========================   =============== ================================================================

birch
^^^^^

For detailed description look at :class:`~sklearn.cluster.Birch` documentation. Following table summarized options available in `Valve` and is a copy of original documentation.

.. tabularcolumns:: |p{2.5cm}|p{2.5cm}|p{9.6cm}|

=========================   =============== ================================================================
Option                      Value type      Description
=========================   =============== ================================================================
threshold                   float           The radius of the subcluster obtained by merging a new sample
                                            and the closest subcluster should be smaller than the threshold.
                                            Otherwise a new subcluster is started.
branching_factor            int             Maximum number of CF subclusters in each node.
n_clusters                  int             Number of clusters after the final clustering step, which
                                            treats the subclusters from the leaves as new samples. By
                                            default, this final clustering step is not performed and the
                                            subclusters are returned as they are.
=========================   =============== ================================================================

kmeans
^^^^^^

For detailed description look at :class:`~sklearn.cluster.KMeans` documentation. The following table summarized options available in `Valve` and is a copy of original documentation.

.. tabularcolumns:: |p{2.5cm}|p{2.5cm}|p{9.6cm}|

=========================   =============== ================================================================
Option                      Value type      Description
=========================   =============== ================================================================
n_clusters                  int             The number of clusters to form as well as the number of
                                            centroids to generate.
max_iter                    int             Maximum number of iterations of the k-means algorithm for a
                                            single run.
n_init                      int             Number of times the k-means algorithm will be run with different
                                            centroid seeds. The final results will be the best output of
                                            n_init consecutive runs in terms of inertia.
init                        str             Method for initialization, defaults to ``k-means++``. Can be
                                            one of following: ``k-means++`` or ``random``.
tol                         float           Relative tolerance with regards to inertia to declare
                                            convergence.
=========================   =============== ================================================================

.. _smoothing_options:

Smooth section
--------------

Section **smooth** supports the following options:

.. tabularcolumns:: |p{2.0cm}|p{2.0cm}|p{10.6cm}|

=========================   =============== ================================================================
Option                      Value type      Description
=========================   =============== ================================================================
method                      str             Smoothing method. Can be one of the following:

                                            * ``window``, (see :class:`~aquaduct.geom.smooth.WindowSmooth`)
                                            * ``mss``, (see :class:`~aquaduct.geom.smooth.MaxStepSmooth`)
                                            * ``window_mss``, (see :class:`~aquaduct.geom.smooth.WindowOverMaxStepSmooth`)
                                            * ``awin``, (see :class:`~aquaduct.geom.smooth.ActiveWindowSmooth`)
                                            * ``awin_mss``, (see :class:`~aquaduct.geom.smooth.ActiveWindowOverMaxStepSmooth`)
                                            * ``dwin``, (see :class:`~aquaduct.geom.smooth.DistanceWindowSmooth`)
                                            * ``dwin_mss``, (see :class:`~aquaduct.geom.smooth.DistanceWindowOverMaxStepSmooth`)
                                            * ``savgol``. (see :class:`~aquaduct.geom.smooth.SavgolSmooth`)
recursive                   int             Number of recursive runs of smoothing method.
window                      int or float    In window-based method defines window size. In plain ``window``
                                            it has to be an int number. In ``savgol`` it has to be odd integer.
step                        int             In step based method defines size of the step.
function                    str             In window based methods defines averaging function. Can be
                                            ``mean`` or ``median``.
polyorder                   int             In ``savgol`` is polynomial order.
=========================   =============== ================================================================
