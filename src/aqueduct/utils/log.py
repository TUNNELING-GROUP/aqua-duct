'''
Created on Dec 15, 2015

@author: tljm
'''

__all__ = ['debug','info','warning','error','critical','pbar']

import logging


from sys import stderr

'''
def pbar(maxval=100):
    # returns new progress bar
    widgets = [pb.Percentage(), ' ',
               pb.Bar(), ' ',
               pb.ETA()]
    return pb.ProgressBar(widgets=widgets,
                          maxval=maxval).start()
'''

#TODO: try to add tqdm as an alternative for progressbarr which does not work very well in ipython


class pbar(object):

    def __init__(self,maxval=100,kind='tqdm'):
        self.__kind = kind

        self.__maxval = maxval

        self.__curval = 0

        if self.__kind == "progressbar":
            import progressbar as pb
            widgets = [pb.Percentage(), ' ',
               pb.Bar(), ' ',
               pb.ETA()]
            self.__pbar = pb.ProgressBar(widgets=widgets,maxval=maxval).start()
        elif self.__kind == "tqdm":
            # it does not work as expected
            from tqdm import tqdm
            self.__pbar = tqdm(leave=True,ascii=True,total=maxval,dynamic_ncols=True,ncols=20)
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
        if self.__kind == "progressbar":
            self.__pbar.finish()



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