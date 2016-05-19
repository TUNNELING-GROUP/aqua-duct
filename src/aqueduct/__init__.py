
"""
Aqueduct - a collection of tools to trace residues in MD simulation.

"""


def version():
    '''
    Returns :mod:`aqueduct` version number.
    
    :return: 3 element tuple of int numbers
    :rtype: tuple
    '''
    return 0,2,7


def version_nice():
    '''
    Returns :mod:`aqueduct` version number as nicely formated string.
    
    :return: string composed on the basis of the number returned by :func:`version`.
    :rtype: str
    '''
    return '.'.join(map(str,version()))


def greetings():
    '''
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
    '''
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
