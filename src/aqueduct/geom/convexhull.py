'''
Created on Dec 4, 2015

@author: tljm
'''

from aqueduct.geom.convexhull_scipy import ConvexHull

#TODO: write proper multiprocessing version of CH
'''
from multiprocessing import Queue,Manager,Lock,Value,Process

# multiprocessing

manager = Manager() 
Q = manager.Queue(4) # queue
'''


def get_queue(no_of_threads):
    
    manager = Manager()
    return manager.Queue(no_of_threads)


class PointWithinConvexHull_Worker(object):

    def __init__(self,chull):
        
        self.chull = chull
        
    def __call__(self,point_queue,results):
        # Provides callable interface. It has to be called with queue of SDF
        # file names and Counter object - here named hits.
         
        # endless loop - waiting for termination signal
        while (True):
            # get current file name out of Queue
            sdf_name = sdf_names_queue.get()
            # None is termination signal
            if sdf_name == None:
                return # terminatio signal detected
            # call match method
            self.match_sdf_file(sdf_name,hits)
            # file is processed - increase number of processed files
            fnr = hits.increment_file()
            # display some useful message
            if hits.display_every_record == 1: # verbose
                message = ("File %d of %d done." + os.linesep) % (fnr,self.nof)
                for nop in xrange(len(self.patt)):
                    message += ("Pattern %d: in total %d matches found so far." + os.linesep) % (nop+1,hits.value(nop=nop))
                logging.info(message[:-1])
            else: # normal output
                logging.info("Found total %d matches so far. File %d of %d done." % (hits.value_total,fnr,self.nof))
