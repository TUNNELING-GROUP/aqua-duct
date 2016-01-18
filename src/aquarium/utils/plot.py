
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

from aquarium.geom import traces


def showit(gen):
    def patched(*args, **kwargs):
        obj = gen(*args, **kwargs)
        plt.show()
        return obj
    return patched

def get_ax3d(fig,sub=111):
    return fig.add_subplot(sub, projection='3d')



class GenericTracePlotter(object):

    @showit
    def init_ax(self):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection='3d') 
        
        self.ax.set_axis_bgcolor('none')
        self.ax.axis('off')
        
        self.fig.subplots_adjust(left=0,bottom=0,right=1,top=1)
        self.fig.set_facecolor('w')

    @showit
    def single_trace(self,coords,color='r',**kwargs):
        coords = np.array(coords)
        self.ax.plot3D(coords[:,0],
                       coords[:,1],
                       coords[:,2],
                       c=color,**kwargs)
    @showit
    def path_trace(self,path,color=('r','g','b'),**kwargs):
        
        for nr,trace in enumerate(path):
            # mid points!
            if len(trace) > 0:
                if nr == 0:
                    if len(path[1]) > 0:
                        midp = np.mean(np.vstack((trace[-1],path[1][0])),0)
                        trace = np.vstack((trace,midp))
                if nr == 1:
                    if len(path[0]) > 0:
                        midp = np.mean(np.vstack((trace[0],path[0][-1])),0)
                        trace = np.vstack((midp,trace))
                    if len(path[2]) > 0:
                        midp = np.mean(np.vstack((trace[-1],path[2][0])),0)
                        trace = np.vstack((trace,midp))
                if nr == 2:
                    if len(path[1]) > 0:
                        midp = np.mean(np.vstack((trace[0],path[1][-1])),0)
                        trace = np.vstack((midp,trace))
                self.single_trace(trace, color=color[nr], **kwargs)
                    
                    
class SimpleProteinPlotter(GenericTracePlotter):

    @showit
    def protein_trace(self,protein,smooth=None,color='c',**kwargs):
        # assumes protein is reader object
        #TODO: iterate over chains?
        bb = protein.parse_selection("protein and backbone")
        coords = bb.get_positions()
        cdiff = traces.diff(coords)

