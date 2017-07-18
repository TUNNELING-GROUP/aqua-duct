# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2016-2017  Tomasz Magdziarz, Alicja Płuciennik, Michał Stolarczyk <info@aquaduct.pl>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Aqua-Duct - a collection of tools to trace residues in MD simulation.
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
    Returns :mod:`aquaduct` version number.

    :return: 3 element tuple of int numbers
    :rtype: tuple
    """
    return 0, 3, 7
    # MDA == 0.15!
    # pbar bug fix
    # window frames support
    # no lock in clui

def version_nice():
    """
    Returns :mod:`aquaduct` version number as nicely formatted string.

    :return: string composed on the basis of the number returned by :func:`version`.
    :rtype: str
    """
    return '.'.join(map(str, version()))


__version__ = version_nice()
__mail__ = 'info@aquaduct.pl'

__author__ = u'Tomasz Magdziarz, Alicja Płuciennik, Michał Stolarczyk'
__author_doc__ = u'Tomasz Magdziarz, Karolina Mitusińska, Sandra Gołdowska, Alicja Płuciennik, Michał Stolarczyk, Magdalena Ługowska, Artur Góra'

def greetings():
    """
    Returns fancy greetings of :mod:`aquaduct`. It has a form of ASCII-like
    graphic. Currently it returns following string::

        ------------------------------------------------
                  ~ ~ ~ A Q U A - D U C T ~ ~ ~
        ################################################
        ####        ########        ########        ####
        ##            ####            ####            ##
        #              ##              ##              #
        #              ##              ##              #
        #              ##              ##              #
        #              ##              ##              #
        ------------------------------------------------

    :return: :mod:`aquaduct` fancy greetings.
    :rtype: str
    """
    greet = '''------------------------------------------------
          ~ ~ ~ A Q U A - D U C T ~ ~ ~
################################################
####        ########        ########        ####
##            ####            ####            ##
#              ##              ##              #
#              ##              ##              #
#              ##              ##              #
#              ##              ##              #
------------------------------------------------'''

    return greet
