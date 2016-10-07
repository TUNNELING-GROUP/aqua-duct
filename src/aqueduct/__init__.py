# -*- coding: utf-8 -*-

"""
Aqueduct - a collection of tools to trace residues in MD simulation.

"""

import logging

try:  # Python 2.7+
    from logging import NullHandler
except ImportError:
    class NullHandler(logging.Handler):
        def emit(self, record):
            pass

logger_name = __name__
logger = logging.getLogger(logger_name)
logger.addHandler(NullHandler())
logger.setLevel(logging.DEBUG)


def version():
    """
    Returns :mod:`aqueduct` version number.

    :return: 3 element tuple of int numbers
    :rtype: tuple
    """
    return 0, 2, 21
    # force-color in pymol vis scripts


def version_nice():
    """
    Returns :mod:`aqueduct` version number as nicely formatted string.

    :return: string composed on the basis of the number returned by :func:`version`.
    :rtype: str
    """
    return '.'.join(map(str, version()))


__version__ = version_nice()
__mail__ = 'info@aquaduct.pl'


def greetings():
    """
    Returns fancy greetings of :mod:`aqueduct`. It has a form of ASCII-like
    graphic. Currently it returns following string::

        ------------------------------------------------
                   ~ ~ ~ A Q U E D U C T ~ ~ ~
        ################################################
        ####        ########        ########        ####
        ##            ####            ####            ##
        #              ##              ##              #
        #              ##              ##              #
        #              ##              ##              #
        #              ##              ##              #
        ------------------------------------------------

    :return: :mod:`aqueduct` fancy greetings.
    :rtype: str
    """
    greet = '''------------------------------------------------
           ~ ~ ~ A Q U E D U C T ~ ~ ~
################################################
####        ########        ########        ####
##            ####            ####            ##
#              ##              ##              #
#              ##              ##              #
#              ##              ##              #
#              ##              ##              #
------------------------------------------------'''

    return greet
