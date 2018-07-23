
*Valve* data
============

Description of data saved by *Valve* in stages 1-4 of calculations and the format of data.

Types of data
-------------

Each of 6 stages of *Valve* calculations can save certain type of data. Each type of data has its unique *name*.
Different stages can save the same type of data by using the same *name*.

Following sections give brief description of saved data.

Stage I **traceable_residues**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Types of data saved in the Stage I:

* `center_of_system` - Coordinates of the center of all atoms of the scope area.

    .. note::

        *Scope* can be redefined in stage II but in the current version
        the center of the system is calculated only in stage I.

* `all_res` - List of IDs of all traced residues in all layers.
* `number_frame_rid_in_object`  - Lists of IDs of traced residues that are in the object in each frame in all layers.


Stage II **raw_paths**
^^^^^^^^^^^^^^^^^^^^^^

Types of data saved in the Stage II:

* `all_res` - List of IDs of all traced residues in all layers.
* `paths` - L found for each traced residues.


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

* `ctypes` - List of cluster-cluster types detected.
* `master_paths` - List of master paths for `ctypes`.
* `master_paths_smooth` - List of smooth master paths for `ctypes`.
* `inls` - List of inlets and clusters.

Common methods for concise data storing
---------------------------------------

Lists of monotically increasing integers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Lists of such types can be represented as collections of continous ranges.
Therefore, they can be stored as list of 2 element tuples (E,T).
Each tuple encodes one continous range. First element *E* is the starting value,
second element *T* is the lenght of the range.

For example, following list of lenght 10::

    [3, 4, 5, 7, 8, 9, 20, 30, 31, 32]

can be stored as matrix 4x2 (8 elements)::

    [[3, 3],
     [7, 3],
     [20, 1],
     [30, 3]]


.. _list_list_monoincr:

Lists of lists of monotically increasing integers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to store lists of the above described lists of monotically increasing integers
in one matrix Nx2 additional array is needed to store sizes of inner lists.
To make it even more efficient, additional array can store number of pairs used to encode
lists instead of true lists lenghts.

For example, following 3 lists, 20 elements in total::

    [[3, 4, 5, 7, 8, 9, 20, 30, 31, 32],
     [3, 4, 5, 6, 7],
     [10, 11, 12, 20, 21]]

can be stored as 7x2 matrix::

    [[3, 3],
     [7, 3],
     [20, 1],
     [30, 3],
     [3, 5],
     [10, 3],
     [20, 2]]

and as additional array of lenght 3::

    [4, 1, 2]

This makes in total 15 elements.

.. _list_arrays_matrices:

Lists of arrays or matrices
^^^^^^^^^^^^^^^^^^^^^^^^^^^

In similar manner lists of arrays (1 dimension) or matrices (having the same all dimenstions but one).

Arrays and matrices can be concatenated. Additional array with lenghts of
concatenated arrays/matrices is needed to allow decoding concatenated data.

For example, 3 matrices Mx3::

    [[1.1, 2.5, 3.7],
     [0.8, 2.3, 3.6],
     [0.3, 2.9, 3.8]]

    [[7.8, 8.5, 4.4],
     [8.0, 9.2, 4.6]]

    [[6.3, 9.5, 3.1],
     [5.9, 9.3, 2.9],
     [5.9, 8.8, 2.7],
     [6.1, 8.5, 2.6]]

can be concateneted to::

    [[1.1, 2.5, 3.7],
     [0.8, 2.3, 3.6],
     [0.3, 2.9, 3.8],
     [7.8, 8.5, 4.4],
     [8.0, 9.2, 4.6],
     [6.3, 9.5, 3.1],
     [5.9, 9.3, 2.9],
     [5.9, 8.8, 2.7],
     [6.1, 8.5, 2.6]]

and additional array with lengths is created::

    [3, 2, 4]


Format of data
--------------

Data *name* is unambiguously related with data *type*. Each data *type* can be saved as one or more
matrices. Each matrix can be of following types:

* int
* float
* str

If several matrices have to be used to encode particular data *type*, each of them has its own name. Matrix name stems from data *name* plus some suffixes. If only one matrix is used its name is equal to data *name*.

Following data *types* can be currently encoded:

* `center_of_system`
* `all_res`
* `number_frame_rid_in_object`
* `paths`
* `spaths`
* `ctypes`
* `master_paths`
* `master_paths_smooth`
* `inls`

`center_of_system`
^^^^^^^^^^^^^^^^^^

.. tabularcolumns:: |p{4.0cm}|p{1.0cm}|p{1.0cm}|p{8.2cm}|


=================   ======  ======   ===================================================================
Matrix name         Shape   Type     Description
=================   ======  ======   ===================================================================
center_of_system    (3,)    float    Cartesian X,Y,Z coordinates of the center of the scope.
=================   ======  ======   ===================================================================

`all_res`
^^^^^^^^^

.. tabularcolumns:: |p{4.0cm}|p{1.0cm}|p{1.0cm}|p{8.2cm}|

=================   ======  ======  ===================================================================
Matrix name         Shape   Type    Description
=================   ======  ======  ===================================================================
all_res.layers      (L,)    int     List of layers.

                                    * *L* is a number of layers.

all_res.layer.N     (S,)    int     List of residues IDs in layer *N*.

                                    * *N* is a layer id.
                                    * *S* is a number of residues in layer *N*.
=================   ======  ======  ===================================================================

`number_frame_rid_in_object`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note::

    IDs of traced residues that are in the object in each frame in all layers are stored as
    :ref:`list_arrays_matrices`: see **nfrio.layer.N.sizes** and **nfrio.layer.N** matrices.

.. tabularcolumns:: |p{4.0cm}|p{1.0cm}|p{1.0cm}|p{8.2cm}|

====================    ======  ======  ===================================================================
Matrix name             Shape   Type    Description
====================    ======  ======  ===================================================================
nfrio.layers.nr         (1,)    int     Number of layers.
nfrio.layer.N.sizes     (F,)    int     Array of numbers of residues indetified in the object area in
                                        frames.

                                        * *N* is a consecutive layer number.
                                        * *F* Number of frames in layer *N*.

nfrio.layer.N           (Q,)    int     IDs of residues indetified in the object area in frames in layer
                                        *N*. It is storead as one list and have to be divieded in to
                                        chunks corresponding to consecutive frames. Sizes of chunks are
                                        stored in nfrio.layer.N matrix

                                        * *Q* is a number of molecules identified in the object in all
                                          frames in layer *N*; it is a sum of nfrio.layer.N.sizes.
====================    ======  ======  ===================================================================


`paths`
^^^^^^^

.. note::

    Lists of frames in which paths are in the object and scope areas are stored as :ref:`list_list_monoincr`:
    see **paths.layer.N.object.sizes**, **paths.layer.N.object**, **paths.layer.N.scope.sizes**, and **paths.layer.N.scope**
    matrices.

.. tabularcolumns:: |p{5.0cm}|p{1.0cm}|p{1.0cm}|p{7.2cm}|

=============================   ========    ======  ===================================================================
Matrix name                     Shape       Type    Description
=============================   ========    ======  ===================================================================
paths.layers                    (L,)        int     List of layers.

                                                    * *L* is a number of layers.

paths.layer.N.min_max_frames    (1,2)       int     Minimal and maximal frame number possible in layer *N*.

                                                    * *N* is a consecutive layer number.

paths.layer.N.names             (P,3)       str     List of residue names corresponding to paths. Each name can
                                                    have only 3 charcters.

                                                    * *N* is a consecutive layer number.
                                                    * *P* is a number of paths in N layer

paths.layer.N.ids               (P,)        int     List of IDs of residues that correspond to paths.

                                                    * *N* is a consecutive layer number.
                                                    * *P* is a number of paths in N layer

paths.layer.N.object.sizes      (P,)        int     Array of numbers of pairs encoding farmes in wich paths are in
                                                    the object area in layer *N*.

                                                    * *N* is a consecutive layer number.

paths.layer.N.scope.sizes       (P,)        int     Array of numbers of pairs encoding farmes in wich paths are in
                                                    the scope area in layer *N*.

                                                    * *N* is a consecutive layer number.

paths.layer.N.object            (PO,2)      int     List of pairs that encodes frames in the object in paths in
                                                    layer *N*.

                                                    * *PO* is a list of pairs used to represent object paths for
                                                      all paths in *N* layer.

paths.layer.N.scope             (PS,2)      int     List of pairs that encodes frames in the scope in paths in
                                                    layer *N*.

                                                    * *PS* is a list of pairs used to represent scope paths for
                                                      all paths in *N* layer.

=============================   ========    ======  ===================================================================


`spaths`
^^^^^^^^

.. note::

    List of frames in which paths are in the object area are stored as :ref:`list_list_monoincr`:
    see **spaths.layer.N.object.sizes** and **spaths.layer.N.object** matrices.

.. tabularcolumns:: |p{4.0cm}|p{1.0cm}|p{1.0cm}|p{8.2cm}|

=============================   ========    ======  ===================================================================
Matrix name                     Shape       Type    Description
=============================   ========    ======  ===================================================================
spaths.layers                   (L,)        int     List of layers.

                                                    * *L* is a number of layers.

spaths.layer.N.names            (P,3)       str     List of residue names corresponding to paths. Each name can
                                                    have only 3 charcters.

                                                    * *N* is a consecutive layer number.
                                                    * *P* is a number of paths in N layer

spaths.layer.N.ids              (P,)        int     List of IDs of residues that correspond to paths.

                                                    * *N* is a consecutive layer number.
                                                    * *P* is a number of paths in N layer

spaths.layer.N.single           (P,)        int     List of flags indicating if spath is a *single* path. If it is
                                                    flag is set `1` otherwise it is set `0`.

                                                    * *N* is a consecutive layer number.
                                                    * *P* is a number of paths in N layer

spaths.layer.N.frames           (P,5)       int     Table decoding spaths. Following columns are in the table:

                                                    #. Starting frame.
                                                    #. End frame.
                                                    #. Lenght of incoming part.
                                                    #. Lenght of object part.
                                                    #. Lenght of outgoing part.

                                                    Passing paths does not need lenghts of parts but are saved in the same was as single paths.

                                                    * *N* is a consecutive layer number.
                                                    * *P* is a number of paths in N layer

spaths.layer.N.object.sizes     (P,)        int     Array of numbers of pairs encoding farmes in wich paths are in
                                                    the scope area in layer *N*.

                                                    * *N* is a consecutive layer number.

spaths.layer.N.object           (PO,2)      int     List of pairs that encodes frames in the object in paths in
                                                    layer *N*.

                                                    * *PO* is a list of pairs used to represent object paths for
                                                      all paths in *N* layer.

=============================   ========    ======  ===================================================================



