'''
Created on Dec 15, 2015

@author: tljm
'''

import logging

import progressbar as pb

def pbar(maxval=100):
    # returns new progress bar
    widgets = [pb.Percentage(),
               pb.Bar(),
               pb.ETA()]
    return pb.ProgressBar(widgets=widgets,
                          maxval=maxval).start()
    