def version():
    return 0,1,0


def version_nice():
    return '.'.join(map(str,version())) + ' 20160321'


def greetings():
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
