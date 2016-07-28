Valve manual
============

Valve application is a kind of driver tha uses :mod:`aqueduct` module to perform analysis of trajectories of selected residues in MD simulation.

Valve invocation
----------------

Once :mod:`aqueduct` module is installed (see :doc:`aqueduct_install`) properly on the machine Valve is available as valve.py command line tool.

Usage
^^^^^

Basic help of Valve usege can be displayed by following command::

    valve.py --help

It should display following information::

    Valve, Aqueduct driver

    optional arguments:
      -h, --help            show this help message and exit
      --dump-template-config
			    Dumps template config file. Supress all other output
			    or actions. (default: False)
      -t THREADS            Limit Aqueduct calculations to given number of
			    threads. (default: None)
      -c CONFIG_FILE        Config file filename. (default: None)


Configuration file template
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Configuration file used by Valeve is of moderate lenght. It can be easily prepared with a template file that cn be printed by Valve.

In order to print configuration file template on the screen type following comand::

    valve.py --dump-template-config

Configuration file template can also be easily saved in to a file with following command::

    valve.py --dump-template-config > config.txt

Where config.txt is a configuration file template.


Valve calculation run
^^^^^^^^^^^^^^^^^^^^^

Once configuration file is ready Valve calculations can be run with a following simple command::

    valve.py -c config.txt

Valve can run some of calulations in parallel. By default it will use all available CPU cores. This is not always desired - limitation of used CPU cores can be done with -t option which limits number of concurent threads used by Valve. If it equals 1 no paralleism is used.

.. note::

    Specyfing number of threads greater then available CPU cores is generally not optimal.

    However, in order to maximize usage of available CPU power it is reccomended to define it as number of cores + 1. The reason is that Valve uses one thread for the main process and the excess over one for processes for parallel calculations.


How does Valve work
-------------------

Application starts with parsing input options. If ``--help`` or ``--dump-template-config`` options are provided appropriate messages are printed on the screen and Valve quits with signal 0.

If config file is provided Valve parse it quickly and regular calculations starts according to its content. Calculations performed by Valve are done in several stages described in the next sections.

Traceable residues
^^^^^^^^^^^^^^^^^^

The first stage finds all residues that should be traced and appends them to the list of tracable resiudes. It is done in a loop over all frames. In each frame residues of interest are searched and appended to the list but only if they are not alredy present on the list.

The search of the residues is done according to he definitions provided by the user. Two requirements have to be met to append residue to the list:

#. The residue have to be foud according to the *Object* definition.
#. The residues have to be withint the *Scope* of interest.

The Object definition encompasses usually the active site of the protein. The scope of interest defines, on the other hand, the boundaries in which residues are traced and is usually defined as protein.

Since :mod:`aqueduct` in its current version uses `MDAnalysis <http://www.mdanalysis.org/>`_ Python module for reading, parsing and searching of MD trajectory data definitions of Object and Scope have to be given as its *Selection Commands*.

Object definition
"""""""""""""""""

Object definition have to comprise of two elements:

#. It have to define residues to trace.
#. It have to define spatial boundaries of the *Object* site.

For example, proper Object definition could be following::

    (resname WAT) and (sphzone 6.0 (resnum 99 or resnum 147))

It defines ``WAT`` as residues that should be traced and defines spatial constrains of the Object site as spherical zone within 6 Angstroms of the center of masses of residue with number 99 and 147.

Scope definition
""""""""""""""""

Scope can be defined in two ways: as Object but with broader boundaries or with the convex hull of selected  molecular object.

In the first case definition is very similar to Object and it have to follow the same limitation. For example, proper Scope definition could be following::

    resname WAT around 2.0 protein

It consequently have to define ``WAT`` as residues of interest and defines spatial constrains as all ``WAT`` residues that are within 2 Angstroms of the protein.

If the scope is defined as the convex hull of selected molecular object (which is recommended), the definition itself have to comprise of this molecular object only. For example ``protein``. In that case the scope is iterpreted as the interior of the convex hull of atoms from the defeinition. Therefore, tracable residues would be in the scope only if they are within the convex hull of atoms of ``protein``.

Raw paths
^^^^^^^^^

The second stage of calculations uses the list of all traceable residues from the first stage and finds coordinates of center of masses for each residue in each frame. As in the first stage it is done in a loop over all frames. For each resiudue in each frame Valve calculates or cheks two things:

#. Is the resiude in the Scope (this is always calculated according to the Scope definition).
#. Is the residue in the Object. This information is calculated in the first stage and can be reused in the second. However, it is also possible to recalculate this data according to the new Object definition.

For each of the tracable resiudues a special Path object is created. If the residue is in the Scope its center of mass is added to the appropriate Path object together with the information if it is in the object or not.

Separate paths
^^^^^^^^^^^^^^

The third stage uses collection of Path objects to create Separate Path objects. Each Path comprise data for one resiude. It may happen that the resiudue enters and leaves the Scope and the Object many times over the entire MD. Each such an event is considered by Valve as a separate path.

Each separate path comprises of three parts:

#. Incoming - Defined as a path that leads from the point in which residue enters the scope and enters the object for the firs time.
#. Object - Defined as a path thet leads from the point in which residue enters the Object for the first time and leaves it for the last time.
#. Outgoing - Defined as a path that leads from the point in which residue leaves the obejct for the last lime and leaves the sope.

Clusterization of inlets
^^^^^^^^^^^^^^^^^^^^^^^^

Each of the separate paths has begining and end. If either of them are at the boundaries onf the Scope they are considered as *Inlets*, i.e. points that mark where the traceable resiudues can enter or leave the Scope. Clusters of inlets, on the other hand, mark endings of tunnels or ways in the system which was simulated in the MD.

Clusterization of inlets is performed in following steps:

#. Initial clusterization. Depending on the method, some of the inlets might to be arranged to any cluster and are considerd as outliers.
#. [Optional] Outliers detecion. Arragemnt of inlets to clusters is sometimes far from optimal. In this step, ilets that do not fit to cluster are detected and annotated as outliers. This step can be executed in two modes:

    #. Automatic mode. Inlet is considered to be an outlier if its distance from the centroid is greater then mean distance + 4 * standard deviation of all distances within the cluster.
    #. Defined threshold. Inlet is considered to be an outlier if its minimal distance from any other point in the cluster is greater then the treshold.

#. [Optional] Reclusterization of outliers. It may happen that the outliers form acctually clusters but it ws not recognized in initial clusteriation. In this step clusterization is executed for outliers only and found clusters are added to the clusters identified in the first step. Rest of the inlets are marked as outliers.

Analysis
^^^^^^^^

Fifth stage of Valve calculations analyses results calculated in stages 1 to 4. Results of the analysis is displayed on the screen or can be save to text file and comprise of following parts:

* Tile and data stamp.
* [Optional] Dump of configuraion options.
* Basic information on traceable residues and separate paths.
    * Number of traceable residues.
    * Number of separate paths.
* Basic information on inlets.
    * Number of inlets.
    * Number of clusters.
    * Are outliers detected.
* Summary of inlets clusters. Table with 5 columns:
    #. **Nr**: Row number, starting from 0.
    #. **Cluster**: ID of the cluster. Outliers have 0.
    #. **Size**: Size of the cluster.
    #. **INCOMING**: Number of inlets corresponding to separate paths that enter the scope.
    #. **OUTGOING**: Number of inlets corresponding to separate paths that leave the scope.
* Summary of separate paths clusters types. Table with 9 columns.
    #. **Nr**: Row number, starting from 0.
    #. **CType**: Separate path Cluster Type.
    #. **Size**: Number of separate paths belonging to Cluster type.
    #. **Inp**: Average lenght of incoming part of the path. If no incoming part is available it is nan.
    #. **InpStd**: Standard deviation of lenght Inp.
    #. **Obj**: Average lenght of object part of the path. If no incoming part is available it is nan.
    #. **ObjStd**: Standard deviation of lenght Inp.
    #. **Out**: Average lenght of outgoing part of the path. If no incoming part is available it is nan.
    #. **OutStd**: Standard deviation of lenght Inp.
* List of separate paths and their properties. Table with 17 columns.
    #. **Nr**: - Row number, starting from 0.
    #. **ID**: - Separate path ID.
    #. **BeginF**: Number of frame in which the path begins.
    #. **InpF**: Number of frame in which path begins Incoming part.
    #. **ObjF**: Number of frame in which path begins Object part.
    #. **OutF**: Number of frame in which path begins Outgoing part.
    #. **EndF**: Number of frame in which the path ends.
    #. **InpL**: Lenght of Incoming part. If no incoming part nan is given.
    #. **ObjL**: Lenght of Object part.
    #. **OutL**: Lenght of Outgoing part. If no outgoing part nan is given.
    #. **InpS**: Average step of Incoming part. If no incoming part nan is given.
    #. **InpStdS**: Standard deviation of InpS.
    #. **ObjS**: Average step of Object part.
    #. **ObjStdS**: Standard deviation of ObjS.
    #. **OutS**: Average step of Outgoing part. If no outgoing part nan is given.
    #. **OutStdS**: Standard deviation of OutS.
    #. **CType**: Cluster type of separate path.

Separate path ID
""""""""""""""""

Separate Paths IDs are composed of two numbers separated by colon. First number is the residue number. Second number is consecutive number of the separate path made by the resiude. Numeration starts with 0.

Cluster Type of separate path
"""""""""""""""""""""""""""""

Each separate paths has two ends: begining and end. Both of them either belong to one of the inlets clusters, or are among ouliers, or are inside the scope. If an end belongs to one of the clusters (including outliers) it has ID of the cluster. If it is inside the scope it has special ID of ``N``. Cluster type is an ID composed of IDs of both ends of separate path separted by colon charcter.

Visualisation
^^^^^^^^^^^^^

Sixth stage of Valve calculations visualizes results calculated in stages 1 to 4. Visualization is done with PyMOL. Valve ustomaticaly starts PyMOL and loads visualisations in to it.
Molecule is loaded as PDB file. Other objects like Inlets clusters or paths are loaded as CGO objects.

Following is a list of objects created in PyMOL (all of them are optional). PyMOL object names given in **bold** text or short explanation is given.

* Selected frame of the simulated system. Object name: *molecule*.
* Inlets clusters, each cluster is a separate object. Object name: **cluster_** followed by cluster annotation: otliers are annotated as Out; regular clusters by ID.
* List of cluster types, raw paths. Each cluster type is a separate object. Object name composed of cluster type (colon replaced by underline) plus **_raw**.
* List of cluster types, smooth paths. Each cluster type is a separate object. Object name composed of cluster type (colon replaced by underline) plus **_smooth**.
* All raw paths. They can be displayed as one object or separated in to Incoming, Object and Outgoing part. Object name: **all_raw**, or **all_raw_in**, **all_raw_obj**, and **all_raw_out**.
* All raw paths inlets arrows. Object name: **all_raw_paths_io**.
* All smooth paths. They can be displayed as one object or separated in to Incoming, Object and Outgoing part. Object name: **all_smooth**, or **all_smooth_in**, **all_smooth_obj**, and **all_smooth_out**.
* All raw paths inlets arrows. Object name: **all_raw_paths_io**.
* Raw paths displayed as separate objects or as one object with several states. Object name: **raw_paths_** plus number of path or **raw_paths** if displayed as one object.
* Smooth paths displayed as separate objects or as one object with several states. Object name: **smooth_paths_** plus number of path or **smooth_paths** if displayed as one object.
* Raw paths arrows displayed as separate objects or as one object with several states. Object name: **raw_paths_io_** plus number of path or **raw_paths_io** if displayed as one object.
* Smooth paths arrows displayed as separate objects or as one object with several states. Object name: **smooth_paths_io_** plus number of path or **smooth_paths_io** if displayed as one object.

Color schemes
"""""""""""""

Inlets clusters are collored automaticaly. Outlaiers are grey.

Incoming parts of paths are red, Outgoing parts are blue. Object parts in case of smooth paths are green and in case of raw paths are green if residue is precisely in the object area or yellow if is leaved object area but it is not in the Outgoing part yet.

Arrows are colored in accordance to paths colors.

Configuration file options
--------------------------

Valve Configuration file is a simple and plain text file. It is similar to INI files comonly used in one of the popular operating systems and is compilant with Python module :mod:`ConfigParser`.

Configuration file comprises of several sections. They can be grupped in to three categories. Names of sections given in **bold** text.

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
^^^^^^^^^^^^^^^^^^

This section allows settings of trajectory data and progress bar type.

Available options
"""""""""""""""""

* ``top`` - Path to Amber topology file.
* ``nc`` - Path to Amber NetCDF file.
* ``pbar`` - Progres bar type. Posible values:
    * ``simple`` - [Default, Recomended] Build in progres bar.

Example
"""""""

::

    [global]
    top = path/to/topology/file.prmtop
    nc = path/to/trajectory/file.nc
    pbar = simple

Common settings for stage sections
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Option **execute**
""""""""""""""""""

All stage sections have ``execute`` option. There are two possible values:

* ``run``
* ``skip``

If ``execute`` is set to ``run`` particular stage is executed, otherwise, it is skipped.

Option **save**
"""""""""""""""

This options allows so save a dump of calculated data on the disk. If ``execute`` is set to ``run`` and ``save`` if a path to the file data calculated in the stage is dumped to the file. If ``save`` is set to ``None`` no data is saved.

Option **load**
"""""""""""""""

If ``execute`` is set to ``skip`` and ``load`` points to the file with saved calculations, saved data is loaded.











