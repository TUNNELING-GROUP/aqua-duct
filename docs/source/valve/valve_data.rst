*Valve* data
============

Description of data saved by *Valve* in all stages of calculations and the format of data.

Types of data
-------------

Each of 6 stages of *Valve* calculations can save certain type of data. Each type of data has its unique *name*.
Different stages can save the same type of data by using the same *name*.

Following sections give description of saved data.

Stage I **traceable_residues**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Types of data saved in the Stage I:

* `center_of_system` - Coordinates of the center of all atoms of the scope area.

    .. note::

        *Scope* can be redefined in stage II but in the current version
        the center of the system is calculated only in stage I.

* `all_res` - List of all traced residues in all layers.
* `number_frame_rid_in_object` IDs of traced residues that are in the object in each frame in all layers.


Stage II **raw_paths**
^^^^^^^^^^^^^^^^^^^^^^

Types of data saved in the Stage II:

* `all_res` - List of all traced residues in all layers.
* `paths` - List of paths found for each traced residues.


Stage III **separate_paths**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Types of data saved in the Stage III:

* `paths` - List of paths found for each traced residues.
* `spaths` - List of separate paths foud for list of paths. Sepeare paths can be of two types:
    * `single` - Paths that overlap with the object area.
    * `passing` - Paths that do not overlap with the obejct area.

Stage IV **inlets_clusterization**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Types of data saved in the Stage IV:

* `ctypes`
* `master_paths`
* `master_paths_smooth`
* `inls`

Stage V **analysis**
^^^^^^^^^^^^^^^^^^^^

Stage VI **visualize**
^^^^^^^^^^^^^^^^^^^^^^

Formats of data
---------------


Method encodes objects into collection of :class:`numpy.ndarray`.
        Objects are identified by name; method does not check if provided value
        is instance of particular classes.

        Following objects are supported:

        * `center_of_system`
        * `all_res`
        * `number_frame_rid_in_object`
        * `paths`
        * `spaths`
        * `ctypes`
        * `master_paths`
        * `master_paths_smooth`
        * `inls`

        Actual structure of the above listed objects is not important and can
        be subject of change in the future. The primary objective of this method
        is to encode *information* stored in these object in a strictly defined manner.

        :param str name: Object name.
        :param value: Acctual object.

        :return: Two element tuple of variable name and :class:`numpy.ndarray`.
        :rtupe: generator
