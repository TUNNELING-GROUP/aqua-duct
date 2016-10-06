Configuration file options
==========================

Valve Configuration file is a simple and plain text file. It has similar structure as INI files commonly used in one of the popular operating systems and is compliant with Python module :mod:`ConfigParser`.

Configuration file comprises of several *sections*. They can be grouped in to three categories. Names of sections given in **bold** text.

#. Global settings:
    * **global**
#. Stages options:
    #. **traceable_residues**
    #. **raw_paths**
    #. **separate_paths**
    #. **inlets_clusterization**
    #. **analysis**
    #. **visualize**
#. Methods options:
    * **smooth**
    * **clusterization**
    * **reclusteriation**

Section **global**
------------------

This section allows settings of trajectory data and some other future global options.

======  =============   ==========================================================================
Option  Default value   Description
======  =============   ==========================================================================
top     None            Path to topology file. Aqueduct supports PDB, PRMTOP, PFS topology files.
trj     None            Path to trajectory file. Aqueduct supports NC and DCD trajectory files.
======  =============   ==========================================================================

.. note::

    Options **top** and **trj** are mandatory.


Common settings of stage sections
---------------------------------

Stages 1-4 which perform calsulations have some common options allowig for execution control and saving/loading data.

========    =================   ===================================================================
Option      Default value       Description
========    =================   ===================================================================
execute     runonce             Option controls stage execution. It can have one of three possible
                                values: ``run``, ``runonce``, and ``skip``. If it is set to ``run``
                                calculations are always performed and if **dump** is set dump file
                                is saved. If it is set to ``runonce`` calculations are performed
                                if there is no dump file specified by **dump** option. If it is
                                present calculations are skiped and data is loaded from the file.
                                If it is set to ``skip`` calculations are skip and if **dump**
                                is set data is loaded from the file.
dump        [dump file name]    File name of dump data. It is used to save results of calculations
                                or to load previously calculated data - this depends on **execute**
                                option. Default value of this option depends on the stage and for
                                stages 1 to 4 is one of the following (listed in order):

                                * 1_traceable_residues_data.dump
                                * 2_raw_paths_data.dump
                                * 3_separate_paths_data.dump
                                * 4_inlets_clusterization_data.dump
========    =================   ===================================================================

Stages 5-6 also uses **execute** option, however, since they do not perform calculations `per se` in stead of **dump** option they use **save**.

========    =================   ===================================================================
Option      Default value       Description
========    =================   ===================================================================
execute     run                 Option controls stage execution. It can have one of three possible
                                values: ``run``, ``runonce``, and ``skip``. If it is set to ``run``
                                or ``runonce`` stage is executed and results is saved according to
                                **save** option. If it is set to ``skip`` stage is skipped.
save        [save file name]    File name for saving results. Default value of this option depends
                                on the stage and for stages 1 to 4 is one of the following
                                (listed in order):

                                * 5_analysis_results.txt
                                * 6_visualize_results.py

                                Stage 6 can save results in two file types:

                                #. As Python script - extension ``.py`` plus companion archive
                                   ``.tar.gz``,
                                #. As PyMOL session - extension ``.pse``.

========    =================   ===================================================================


Stage **traceable_residues**
----------------------------

=================   ==============  ================================================================
Option              Default value   Description
=================   ==============  ================================================================
scope               None            Definition of *Scope* of interest. See also
                                    :ref:`scope_definition`.
scope_convexhull    True            Flag to set if the *Scope* is direct or convex hull definition.
object              None            Definition of *Object* of interest. See also
                                    :ref:`object_definition`.
=================   ==============  ================================================================


.. note::

    Options **scope** and **object** are mandatory.


Stage **raw_paths**
-------------------

This stage also requires definition of the *Scope* and *Object*. If appropriate settings are not given, settings from the previous stage are used.

=====================   ==============  ================================================================
Option                  Default value   Description
=====================   ==============  ================================================================
scope                   None            Definition of *Scope* of interest. See also
                                        :ref:`scope_definition`. If ``None`` value form previous stage
                                        is used.
scope_convexhull        None            Flag to set if the *Scope* is direct or convex hull definition.
                                        If ``None`` value form previous stage is used.
object                  None            Definition of *Object* of interest. See also
                                        :ref:`object_definition`. If ``None``, value form the previous
                                        stage is used
clear_in_object_info    False           If it is set to ``True`` information on occupation of *Object*
                                        site by traceable residues calculated in the previous stage is
                                        cleared and have to be recalculated. This is useful if
                                        definition of *Object* was changed.
=====================   ==============  ================================================================


Stage **separate_paths**
------------------------

=====================   ==============  ================================================================
Option                  Default value   Description
=====================   ==============  ================================================================
discard_empty_paths     True            If set to ``True`` empty paths are discarded.
sort_by_id              True            If set to ``True`` separate paths are sorted by ID. Otherwise
                                        they are sorted in order of apparance.
apply_smoothing         False           If set to ``True`` smooth paths are precalculated according to
                                        **smooth** setting. This speed up access to smooth paths in
                                        later stages but makes dump data much bigger.
apply_soft_smoothing    True            If set to ``True`` raw paths are replaced by smooth paths
                                        calculated according to **smooth** section.
discard_short_paths     1               This option allows to discard paths that are shorter then the
                                        threshold.
auto_barber             None            This option allows to select molecular entity used in Auto
                                        Barber procedure. See also :ref:`auto_barber_procedure` and
                                        :meth:`~aqueduct.traj.paths.GenericPaths.barber_with_spheres`.
=====================   ==============  ================================================================


Stage **inlets_clusterization**
-------------------------------

=====================   ==============  ================================================================
Option                  Default value   Description
=====================   ==============  ================================================================
recluster_outliers      False           If set to ``True`` reclusterization of outliers is executed
                                        according to the method defined in **reclusterization** section.
detect_outliers         False           If set detection of outliers is executed. It could be set as a
                                        floating point distance threshold or set tu ``Auto``. See
                                        :ref:`clusterization_of_inlets` for more details.
singletons_outliers     False           Maximal size of cluster to be considered as outliers. If set to
                                        number > 0 clusters of that size are removed and their objects
                                        are moved to outliers. See :ref:`clusterization_of_inlets` for
                                        more details.
max_level               5               Maximal number of recursive clusterization levels.
create_master_paths     False           If set to ``True`` master paths are created (fast CPU and big
                                        RAM recommended; 50k frames long simulation may need ca 20GB of
                                        memory)
=====================   ==============  ================================================================

Stage **analysis**
------------------

=====================   ==============  ================================================================
Option                  Default value   Description
=====================   ==============  ================================================================
dump_config             True            If set to ``True`` configuration options, as seen by Valve, are
                                        added to the head of results.
=====================   ==============  ================================================================


Stage **visualize**
-------------------

=====================   ================    ========================================================================================
Option                  Default value       Description
=====================   ================    ========================================================================================
simply_smooths          RecursiveVector     Option indicates linear simplification method to be used in
                                            plotting smooth paths. Simplification removes points which do
                                            not (or almost do not) change the shape of smooth path.
                                            Possible choices are:

                                            * ``RecursiveVector`` (see :class:`~aqueduct.geom.traces.LinearizeRecursiveVector`),
                                            * ``HobbitVector`` (see :class:`~aqueduct.geom.traces.LinearizeHobbitVector`),
                                            * ``OneWayVector`` (see :class:`~aqueduct.geom.traces.LinearizeOneWayVector`),
                                            * ``RecursiveTriangle`` (see :class:`~aqueduct.geom.traces.LinearizeRecursiveTriangle`),
                                            * ``HobbitTriangle`` (see :class:`~aqueduct.geom.traces.LinearizeHobbitTriangle`),
                                            * ``OneWayTriangle`` (see :class:`~aqueduct.geom.traces.LinearizeOneWayTriangle`).

                                            Optionally name of the method can be followed by a threshold
                                            value in parentheses, ie ``RecursiveVector(0.05)``. For sane
                                            values of thresholds see appropriate documentation of each method.
                                            Default values works well. This option is not case sensitive.
                                            It is recommended to use default method or ``HobbitVector`` method.
all_paths_raw           False               If True produces one object in PyMOL that holds all paths
                                            visualized by raw coordinates.
all_paths_smooth        False               If True produces one object in PyMOL that holds all paths
                                            visualized by smooth coordinates.
all_paths_split         False               If is set True objects produced by **all_paths_raw** and
                                            **all_paths_smooth** are split into Incoming, Object, and
                                            Outgoing parts and visualized as three different objects.
all_paths_raw_io        False               If set True arrows pointing beginning and end of paths are
                                            displayed oriented accordingly to raw paths orientation.
all_paths_smooth_io     False               If set True arrows pointing beginning and end of paths are
                                            displayed oriented accordingly to smooth paths orientation.
paths_raw               False               If set True raw paths are displayed as separate objects or as
                                            one object with states corresponding to number of path.
paths_smooth            False               If set True smooth paths are displayed as separate objects or
                                            as one object with states corresponding to number of path.
paths_raw_io            False               If set True arrows indicating beginning and and of paths,
                                            oriented accordingly to raw paths, are displayed as separate
                                            objects or as one object with states corresponding to number
                                            of paths.
paths_smooth_io         False               If set True arrows indicating beginning and and of paths,
                                            oriented accordingly to smooth paths, are displayed as separate
                                            objects or as one object with states corresponding to number
                                            of paths.
paths_states            False               If True objects displayed by **paths_raw**, **paths_smooth**,
                                            **paths_raw_io**, and **paths_smooth_io** are displayed as one
                                            object with with states corresponding to number of paths.
                                            Otherwise they are displayed as separate objects.
ctypes_raw              False               Displays raw paths in a similar manner as non split
                                            **all_paths_raw** but each cluster type is displayed in
                                            separate object.
ctypes_smooth           False               Displays smooth paths in a similar manner as non split
                                            **all_paths_smooth** but each cluster type is displayed in
                                            separate object.
show_molecule           False               If is set to selection of some molecular object in the system,
                                            for example to ``protein``, this object is displayed.
show_molecule_frames    0                   Allows to indicate which frames of object defined by
                                            **show_molecule** should be displayed. It is possible to set
                                            several frames. In that case frames would be displayed as
                                            states.
show_chull              False               If is set to selection of some molecular object in the system,
                                            for example to ``protein``, convex hull of this object is
                                            displayed.
show_chull_frames       0                   Allows to indicate for which frames of object defined by
                                            **show_chull** convex hull should be displayed. It is possible
                                            to set several frames. In that case frames would be displayed
                                            as states.
=====================   ================    ========================================================================================


.. note::

    Possibly due to limitations of :mod:`MDAnalysis` only whole molecules can be displayed. If **show_molecule** is set to ``backbone`` complete protein will be displayed any way. This may change in future version of :mod:`MDAnalysis` and or :mod:`aqueduct`.

.. note::

    If several frames are selected they are displayed as states which may interfere with other PyMOL objects displayed with several states.

.. note::

    If several states are displayed protein tertiary structure data might be lost. This seems to be limitation of either :mod:`MDAnalysis` or PyMOL.

.. _clusterization_options:

Clusterization sections
-----------------------

Default section for definition of clusterization method is named **clusterization** and default section for reclusterization method definition is named **reclusterization**. All clusterization sections shares some common options. Other options depends on the method.

=========================   =============== ================================================================
Option                      Default value   Description
=========================   =============== ================================================================
method                      meanshift or    Name of clasteriation method. It have to be one of the
                            dbscan          following: dbscan, affprop, meanshift, birch, kmeans. Default
                                            value depends if it is **clusteriation** section (meanshift) or
                                            **reclusterization** section (dbscan).
recursive_clusterization    clusterization  If it is set to name of some section that holds clusterization
                            or None         method settings this method will be called in the next
                                            recursion of clusteriation. Default value for
                                            **reclusterization** is None.
recursive_threshold         None            Allows to set threshold of that excludes clusters of certain
                                            size from reclusterization. Value of this option comprises of
                                            `operator` and `value`. Operator can be one of the following:
                                            >, >=, <=, <. Value have to be expressed as floating number and
                                            it have to be in the range of 0 to 1.
=========================   =============== ================================================================

dbscan
^^^^^^

For detailed description look at :class:`sklearn.cluster.DBSCAN` documentation. Following table summarized options available in `Valve` and is a copy of original documentation.

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

                                            * ``cityblock``,
                                            * ``cosine``,
                                            * ``euclidean``,
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

For detailed description look at :class:`~sklearn.cluster.AffinityPropagation` documentation. Following table summarized options available in `Valve` and is a copy of original documentation.

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

=========================   =============== ================================================================
Option                      Value type      Description
=========================   =============== ================================================================
threshold                   float           The radius of the subcluster obtained by merging a new sample
                                            and the closest subcluster should be lesser than the threshold.
                                            Otherwise a new subcluster is started.
branching_factor            int             Maximum number of CF subclusters in each node.
n_clusters                  int             Number of clusters after the final clustering step, which
                                            treats the subclusters from the leaves as new samples. By
                                            default, this final clustering step is not performed and the
                                            subclusters are returned as they are.
=========================   =============== ================================================================

kmeans
^^^^^^

For detailed description look at :class:`~sklearn.cluster.KMeans` documentation. Following table summarized options available in `Valve` and is a copy of original documentation.

=========================   =============== ================================================================
Option                      Value type      Description
=========================   =============== ================================================================
n_clusters                  int             The number of clusters to form as well as the number of
                                            centroids to generate.
max_iter                    int             Maximum number of iterations of the k-means algorithm for a
                                            single run.
n_init                      int             Number of time the k-means algorithm will be run with different
                                            centroid seeds. The final results will be the best output of
                                            n_init consecutive runs in terms of inertia.
init                        str             Method for initialization, defaults to ``k-means++``. Can be
                                            one of following: ``k-means++`` or ``random``.
tol                         float           Relative tolerance with regards to inertia to declare
                                            convergence.
=========================   =============== ================================================================

Smooth section
--------------

Section **smooth** supports following options:

=========================   =============== ================================================================
Option                      Value type      Description
=========================   =============== ================================================================
method                      str             Smoothing method. Can be one of the following:

                                            * ``window``, (see :class:`~aqueduct.geom.smooth.WindowSmooth`)
                                            * ``mss``, (see :class:`~aqueduct.geom.smooth.MaxStepSmooth`)
                                            * ``window_mss``, (see :class:`~aqueduct.geom.smooth.WindowOverMaxStepSmooth`)
                                            * ``awin``, (see :class:`~aqueduct.geom.smooth.ActiveWindowSmooth`)
                                            * ``awin_mss``, (see :class:`~aqueduct.geom.smooth.ActiveWindowOverMaxStepSmooth`)
                                            * ``dwin``, (see :class:`~aqueduct.geom.smooth.DistanceWindowSmooth`)
                                            * ``dwin_mss``. (see :class:`~aqueduct.geom.smooth.DistanceWindowOverMaxStepSmooth`)
recursive                   int             Number of recursive runs of smoothing method.
window                      int or float    In window based method defines window size. In plain ``window``
                                            it has to be int number.
step                        int             In step based method defines size of the step.
function                    str             In window based methods defines averaging function. Can be
                                            ``mean`` or ``median``.
=========================   =============== ================================================================
