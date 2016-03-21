from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import colorConverter
import matplotlib.pyplot as plt
import numpy as np

from aqueduct.geom import traces
# matplotlib specific
# mpl color converter
#cc = lambda c,alpha=1.0 : colorConverter.to_rgba(c,alpha=alpha)
cc = lambda c,alpha=1.0 : colorConverter.to_rgb(c)

def showit(gen):
    def patched(*args, **kwargs):
        obj = gen(*args, **kwargs)
        plt.show()
        return obj
    return patched

def get_ax3d(fig,sub=111):
    return fig.add_subplot(sub, projection='3d')



class SimpleTracePlotter(object):

    def plot_line(self,coords,color,**kwargs):
        raise NotImplementedError('This is base class.')


    '''
    def path_trace(self,path,color=('r','g','b'),
                   plot_in=True,
                   plot_object=True,
                   plot_out=True,
                   **kwargs):
        raise NotImplementedError('This is base class.')
    '''


    def single_trace(self,coords,color='r',**kwargs):
        # coords is a trace
        # color is a single color
        color = cc(color)
        coords = np.array(coords)
        # call plot_line
        self.plot_line(coords,color,**kwargs)


    def path_trace(self,path,color=('r','g','b'),
                   plot_in=True,
                   plot_object=True,
                   plot_out=True,
                   **kwargs):
        # path is a tuple of length 3, its elements represent in,object, out parts of path
        # color is a tuple of length 3, its elements correspond to colors of consecutive path parts
        color = map(cc,color)
        for nr,trace in enumerate(traces.midpoints(path)):
            # mid points!
            if len(trace) > 0:
                if (nr == 0 and plot_in) or (nr == 1 and plot_object) or (nr == 2 and plot_out):
                    self.single_trace(trace, color=color[nr], **kwargs)
                    
                    
class SimpleProteinPlotter(SimpleTracePlotter):

    def protein_trace(self,protein,smooth=None,color=('c','m','y'),**kwargs):
        # assumes protein is reader object
        #TODO: iterate over chains?
        bb = protein.parse_selection("protein and backbone")
        coords = bb.get_positions()
        cdiff = traces.diff(coords)
        split = np.argwhere(cdiff > 2.5) #TODO: magic constant
        ns = len(split)
        if ns == 0:
            self.single_trace(smooth(coords),color=color[0],**kwargs) #TODO: color conversion is buggy
        else:
            split.shape = (ns,)
            for nr,csplit in enumerate([0]+split.tolist()):
                cc = color[nr%len(color)]
                if nr == ns:
                    scoords = coords[csplit+1:]
                elif nr == 0:
                    scoords = coords[csplit:split[nr]]
                else:
                    scoords = coords[csplit+1:split[nr]]
                if smooth:
                    scoords = smooth(scoords)
                self.single_trace(scoords,color=color,**kwargs)


class SimplePathPlotter(SimpleTracePlotter):

    def single_path_traces(self,spaths,smooth=None,color=('r','g','b'),**kwargs):
        for spath in spaths:
            self.path_trace(spath.get_smooth_coords(smooth),color=color,**kwargs)



class MPLTracePlotter(SimplePathPlotter, SimpleProteinPlotter):

    @showit
    def init_ax(self):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection='3d')

        self.fig.subplots_adjust(left=0,bottom=0,right=1,top=1)
        self.fig.set_facecolor('w')

        self.ax.set_axis_bgcolor('none')
        self.ax.axis('off')

    @showit
    def plot_line(self,coords,color,**kwargs):
        self.ax.plot3D(coords[:,0],
                       coords[:,1],
                       coords[:,2],
                       c=color,**kwargs)



