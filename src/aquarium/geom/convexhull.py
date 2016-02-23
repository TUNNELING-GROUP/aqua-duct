'''
Created on Dec 4, 2015

@author: tljm
'''

from aquarium.geom.convexhull_scipy import ConvexHull

#TODO: write proper multiprocessing version of CH
'''
from multiprocessing import Queue,Manager,Lock,Value,Process

# multiprocessing

manager = Manager() 
Q = manager.Queue(4) # queue
'''


'''
def get_queue(no_of_threads):
    
    manager = Manager()
    return manager.Queue(no_of_threads)


class PointWithinConvexHull_Worker(object):
'''