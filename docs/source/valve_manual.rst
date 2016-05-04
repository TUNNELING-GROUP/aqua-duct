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

    Specyfing number of threads greater then available CPU cores is not optimal.

Configuration file options
--------------------------





