'''
Created on Dec 15, 2015

@author: tljm
'''

__all__ = ['debug','info','warning','error','critical','pbar']

import logging
import time

from os import linesep
from sys import stderr


def smart_time_string(s, rl=0):
    output = ''
    rl += 1
    if rl > 2:
        output = ''
    else:
        t = 1.1
        # seconds
        if s < t * 60:
            output = "%2.2d s" % s
        elif s < t * 3600:
            output = ("%2.2d m" % (int(s) / 60)) + ' ' + smart_time_string(int(s) % 60, rl).strip()
        elif s < t * 3600 * 24:
            output = ("%2.2d h" % (int(s) / 3600)) + ' ' + smart_time_string(int(s) % 3600, rl).strip()
        elif True:  # s < t*3600*24*365:
            output =  ("%d d" % (int(s) / (3600 * 24))) + ' ' + smart_time_string(int(s) % (3600 * 24), rl).strip()
    return (output+" "*10)[:10]

class SimpleProgressBar(object):

    def __init__(self,maxval=None):
        assert isinstance(maxval,(int,long))
        assert maxval > 0

        self.maxval = maxval
        self.current = 0

        self.overrun_notice = True
        self.overrun = False
        
        self.begin = time.time()
        self.tcurrent = self.begin
        self.show()

    def ETA(self):
        if self.current == 0:
            return '?'
        diff = self.tcurrent - self.begin
        periteration = diff/self.current
        expected = periteration*self.maxval
        eta = periteration*(self.maxval-self.current)

        return smart_time_string(eta)

    def percent(self):
        percent = float(self.current)/float(self.maxval)*100
        return percent


    def show(self):
        percent = self.percent()
        if percent > 100:
            if self.overrun_notice:
                stderr.write(linesep)
                self.overrun_notice = False
                self.overrun = True
            stderr.write("\r%d iterations out of %d. Total time: %s" % (self.current,self.maxval,self.ttime()))
        elif not self.overrun:
            stderr.write("\r%3d%% ETA: %s" % (self.percent(),self.ETA()))

    def update(self,step):
        if step > 0:
            self.current += 1
        self.tcurrent = time.time()
        self.show()

    def ttime(self):
        return smart_time_string(self.tcurrent-self.begin)

    def finish(self):
        self.update(0)
        stderr.write(linesep)
        stderr.write("Total time: %s" % self.ttime())
        stderr.write(linesep)

#TODO: try to add tqdm as an alternative for progressbarr which does not work very well in ipython


class pbar(object):

    def __init__(self,maxval=100,kind='simple'):
        self.__kind = kind

        self.__maxval = maxval

        self.__curval = 0

        if self.__kind == "simple":
            self.__pbar = SimpleProgressBar(maxval=maxval)
        if self.__kind == "progressbar":
            import progressbar as pb
            widgets = [pb.Percentage(), ' ',
               pb.Bar(), ' ',
               pb.ETA()]
            self.__pbar = pb.ProgressBar(widgets=widgets,maxval=maxval).start()
        elif self.__kind == "tqdm":
            # it does not work as expected
            from tqdm import tqdm
            self.__pbar = tqdm(leave=True,ascii=True,total=maxval)
        elif self.__kind == "pyprind":
            import pyprind
            self.__pbar = pyprind.ProgBar(maxval)

    def update(self,val):
        if self.__kind in ["tqdm", "pyprind"]:
            if val < 1:
                return
        if self.__kind in ["tqdm", "pyprind"]:
            self.__pbar.update()
        else:
            self.__pbar.update(val)

    def finish(self):
        if self.__kind in ["progressbar","simple"]:
            self.__pbar.finish()
        if self.__kind == "tqdm":
            self.__pbar.close()



level = logging.WARNING
#format = linesep+'AQUARIUM:%(levelname)1.1s:[%(module)s|%(funcName)s@s%(lineno)d]:'+linesep+'%(message)s'
format = linesep+'AQUARIUM:%(levelname)s:[%(module)s|%(funcName)s@s%(lineno)d]:'+linesep+'%(message)s'

logging.basicConfig(format=format, level=level)

debug = logging.debug
info = logging.info
warning = logging.warning
error = logging.error
critical = logging.critical

class fbm(object):
    # feedback message
    def __init__(self,info):
        message(info+'...',cont=True)
    def __enter__(self):
        pass
    def __exit__(self, type, value, traceback):
        message("OK.")

def message(mess,cont=False):
    if cont:
        print >> stderr, mess,
    else:
        print >> stderr, mess