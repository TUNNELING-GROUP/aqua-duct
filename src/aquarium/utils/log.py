'''
Created on Dec 15, 2015

@author: tljm
'''

__all__ = ['debug','info','warning','error','critical','pbar']

import logging

import progressbar as pb

from sys import stderr

def pbar(maxval=100):
    # returns new progress bar
    widgets = [pb.Percentage(), ' ',
               pb.Bar(), ' ',
               pb.ETA()]
    return pb.ProgressBar(widgets=widgets,
                          maxval=maxval).start()

#TODO: try to add tqdm as an alternative for progressbarr which does not work very well in ipython


level = logging.CRITICAL
format = 'AQUARIUM:%(levelname)1.1s:[%(module)s|%(funcName)s@s%(lineno)d]: %(message)s'

logging.basicConfig(format=format, level=level)

debug = logging.debug
info = logging.info
warning = logging.warning
error = logging.error
critical = logging.critical

def message(mess):
    print >> stderr, mess