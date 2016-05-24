Configuration file options
==========================

Valve Configuration file is a simple and plain text file. It is similar to INI files commonly used in one of the popular operating systems and is compliant with Python module :mod:`ConfigParser`.

Configuration file comprises of several *sections*. They can be grouped in to three categories. Names of sections given in **bold** text.

#. Global settings:
    * **global**
#. Stages options:
    * **traceable_residues**
    * **raw_paths**
    * **separate_paths**
    * **inlets_clusterization**
    * **analysis**
    * **visualize**
#. Methods options:
    * **smooth**
    * **clusterization**
    * **reclusteriation**

Section **global**
------------------

This section allows settings of trajectory data and progress bar type.

Available options
^^^^^^^^^^^^^^^^^

* ``top`` - Path to Amber topology file.
* ``nc`` - Path to Amber NetCDF file.
* ``pbar`` - Progres bar type. Possible values:
    * ``simple`` - [Default, Recommended] Build in progress bar.

Example
^^^^^^^

::

    [global]
    top = path/to/topology/file.prmtop
    nc = path/to/trajectory/file.nc
    pbar = simple

Common settings of stage sections
---------------------------------

Option **execute**
^^^^^^^^^^^^^^^^^^

All stage sections have ``execute`` option which decides if the stage is executed or skipped. There are two possible values of ``execute`` option: ``run``, and ``skip``.

If ``execute`` is set to ``run`` the stage is executed and results of calculations can be optionally saved.

If ``execute`` is set to ``skip`` execution of the stage is skipped. Results of calculations saved previously can be optionally loaded.

Option **save**
^^^^^^^^^^^^^^^

This options allows to save a dump of calculated data on the disk. If ``execute`` is set to ``run`` and ``save`` is set to file name results are saved as gziped pickled dump.

In case of **analysis** and **visualize** sections this setting has slightly different function. Results of **analysis** section as saved to the file pointed by **save** as plain text comprising of tables and summaries. Stage **visualize**, on the other hand, uses **save** option to save PyMOL session.

Option **load**
^^^^^^^^^^^^^^^

This options allows to read previously saved dump of results. It is used only if ``ececute`` is set to ``skip`` and is valid only for **traceable_residues**, **raw_paths**, **separate_paths**, and **inlets_clusterization** stages.

Stage **traceable_residues**
----------------------------

Option **scope**
^^^^^^^^^^^^^^^^

Definition of *Scope* of interest. See also :ref:`scope_definition`.

.. note::

    This option is mandatory.

Option **scope_convexhull**
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Flag to set if the *Scope* is direct or convex hull definition.


Option **object**
^^^^^^^^^^^^^^^^^

Definition of *Object* of interest. See also :ref:`object_definition`.

.. note::

    This option is mandatory.


Stage **raw_paths**
-------------------

This stage also requires definition of the *Scope* and *Object*. If appropriate settings are not given, settings from the previous stage are used.

Option **clear_in_object_info**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If it is set to ``True`` information on occupation of *Object* site by traceable residues calculated in the previous stage are cleared and have to be recalculated. This is useful if definition of *Object* is changed.

Stage **separate_paths**
------------------------

Option **discard_empty_paths**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If set to ``True`` empty paths are discarded.

Option **sort_by_id**
^^^^^^^^^^^^^^^^^^^^^

If set to ``True`` separate paths are sorted by ID.


Option **apply_smoothing**
^^^^^^^^^^^^^^^^^^^^^^^^^^

If set to ``True`` smooth paths are precalculated according to **smooth** setting.
This speed up access to smooth paths in later stages but makes dump data much bigger.


Option **apply_soft_smoothing**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If set to ``True`` raw paths are replaced by smooth paths calculated according to **smooth** section.

Option **discard_short_paths**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This option allows to discard paths that are shorter then the threshold.

Stage **inlets_clusterization**
-------------------------------

Option **recluster_outliers**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If set to ``True`` reclusterization of outliers is executed according to the method defined in **reclusterization** section.

Option **detect_outliers**
^^^^^^^^^^^^^^^^^^^^^^^^^^

If set detection of outliers is executed. See :ref:`clusterization_of_inlets` for more details.

Option **singletons_outliers**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Maximal size of cluster to be considered as outliers. If set to number > 0 clusters of that size are removed and their objects are moved to outliers. See :ref:`clusterization_of_inlets` for more details.


Stage **analysis**
------------------

Option **dump_config**
^^^^^^^^^^^^^^^^^^^^^^

If set to ``True`` configuration options, as seen by Valve, are added to the head of results.


Stage **visualize**
-------------------

Option **simply_smooths**
^^^^^^^^^^^^^^^^^^^^^^^^^

If set to float number simplification of smooth paths is applied.
Simplification removes points which do not (or almost do not) change the shape of smooth path. For more details see :ref:`Recursive Vector Linearization <simply_smooths_details>`.

Option **all_paths_raw**
^^^^^^^^^^^^^^^^^^^^^^^^

If True produces one object in PyMOL that holds all paths visualized by raw coordinates.

Option **all_paths_smooth**
^^^^^^^^^^^^^^^^^^^^^^^^^^^

If True produces one object in PyMOL that holds all paths visualized by smooth coordinates.

Option **all_paths_split**
^^^^^^^^^^^^^^^^^^^^^^^^^^

If is set True objects produced by **all_paths_raw** and **all_paths_smooth** are split into Incoming, Object, and Outgoing parts and visualized as three different objects.

Options **all_paths_raw_io** and **all_paths_smooth_io**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If set True arrows pointing beginning and end of paths are displayed oriented accordingly to raw or smooth paths.

Option **paths_raw**
^^^^^^^^^^^^^^^^^^^^

If set True raw paths are displayed as separate objects or as one object with states corresponding to number of path.

Option **paths_raw**
^^^^^^^^^^^^^^^^^^^^

If set True smooth paths are displayed as separate objects or as one object with states corresponding to number of path.

Option **paths_raw_io**
^^^^^^^^^^^^^^^^^^^^^^^

If set True arrows indicating beginning and and of paths, oriented accordingly to raw paths, are displayed as separate objects or as one object with states corresponding to number of paths.

Option **paths_smooth_io**
^^^^^^^^^^^^^^^^^^^^^^^^^^

If set True arrows indicating beginning and and of paths, oriented accordingly to smooth paths, are displayed as separate objects or as one object with states corresponding to number of paths.

Option **paths_states**
^^^^^^^^^^^^^^^^^^^^^^^

If True objects displayed by **paths_raw**, **paths_smooth**, **paths_raw_io**, and **paths_smooth_io** are displayed as one object with with states corresponding to number of paths. Otherwise they are displayed as separate objects.

Option **ctypes_raw**
^^^^^^^^^^^^^^^^^^^^^

Displays raw paths in a similar manner as non split **all_paths_raw** but each cluster type is displayed in separate object.

Option **ctypes_smooth**
^^^^^^^^^^^^^^^^^^^^^^^^

Displays smooth paths in a similar manner as non split **all_paths_smooth** but each cluster type is displayed in separate object.


Option **show_molecule**
^^^^^^^^^^^^^^^^^^^^^^^^

If is set to selection of some molecular object in the system, for example to ``protein``, this object is displayed.

.. note::

    Possibly due to limitations of :mod:`MDAnalysis` only whole molecules can be displayed. If **show_molecule** is set to ``backbone`` complete protein will be displayed any way. This may change in future version of :mod:`MDAnalysis` and or :mod:`aqueduct`.

Option **show_molecule_frames**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Allows to indicate which frames of object defined by **show_molecule** should be displayed. It is possible to set several frames. In that case frames would be displayed as states.

.. note::

    If several frames are selected they are displayed as states which may interfere with other PyMOL objects displayed with several states.

.. note::

    If several states are displayed protein tertiary structure data might be lost. This seems to be limitation of either :mod:`MDAnalysis` or PyMOL.
