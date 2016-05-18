import numpy as np
import scipy.linalg

import matplotlib as mpl
mpl.rcParams["backend"] = "Qt4Agg"
mpl.interactive(True)
from matplotlib import cm
import pylab as pl

class Center(object):
    def __init__(self,X):
        self.mean = np.mean(X,0)
    def __call__(self,X):
        return X-self.mean
        
class Normalize(object):
    def __init__(self,X):
        self.std = np.std(X,0)
    def __call__(self,X):
        return X/self.std


class Standartize(Center,Normalize):
    def __init__(self,X):
        Center.__init__(self,X)
        Normalize.__init__(self,X)
    def __call__(self,X):
       return  Normalize.__call__(self,Center.__call__(self,X))

def dropna(X):
    return X.dropna(1,how='all').dropna()


class PCA(object):

    def __init__(self,X):
        # SVD
        #self.U, self.d, self.Pt = np.linalg.svd(X, full_matrices=False)
        self.U, self.d, self.Pt = scipy.linalg.svd(X, full_matrices=False)
        assert np.all( self.d[:-1] >= self.d[1:] )  # sorted
        # PCA
        self.T = self.U * self.d
        self.eigen = self.d**2
        self.sumvariance = np.cumsum(self.eigen)
        self.sumvariance /= self.sumvariance[-1]

    @property
    def P(self):
        return self.Pt.T


    @staticmethod
    def prepro(X):
        X = dropna(X)
        return dropna(Standartize(X)(X))

     

import statsmodels.api as sm

def get_XY(deskryptory):
    return dropna(deskryptory).drop('Price_1g',axis=1),dropna(deskryptory)['Price_1g']

def get_OLS_model_results(X,Y):
    model = sm.OLS(Y,X)
    return model.fit()
    
def get_consecutive_components(T):
    for pc in range(1,T.shape[1]):
        yield T[:,:pc]

def do_PCR(deskryptory):
    X,Y = get_XY(deskryptory)
    X_pca = PCA(PCA.prepro(X))
    results = []
    for nr,TT in enumerate(get_consecutive_components(X_pca.T)):
        print nr+1,get_OLS_model_results(TT,Y).rsquared
        results.append((nr+1,get_OLS_model_results(TT,Y).rsquared))
    return results



################################################################################
# plots


def set_axis(X,Y,scale=1.2):
    #max_valx = max(X)*scale*scale
    #max_valy = max(Y)*scale*scale
    #min_valx = min(X)*scale*scale
    #min_valy = min(Y)*scale*scale
    
    max_valx = max(max(X)*scale*scale,0)
    max_valy = max(max(Y)*scale*scale,0)
    min_valx = min(min(X)*scale*scale,0)
    min_valy = min(min(Y)*scale*scale,0)

    max_val = max(max_valx,max_valy)
    min_val = min(min_valx,min_valy)
    
    pl.axis((min_valx/scale,max_valx/scale,min_valy/scale,max_valy/scale))
    #pl.axis((min_val/scale,max_val/scale,min_val/scale,max_val/scale))


def plot_PCA(D, names=None, extreme=True, PC1=1, PC2=2, rays=False, origin=False, format=".", scale=0.2, **args):
    # D - data
    # N - names
    X = D[:,PC1-1]
    Y = D[:,PC2-1]
    # rays
    if rays == True: # origin
        for x,y in zip(X,Y):
            pl.plot([0,x],[0,y],'-',color=[0.85,0.9,0.85])
    # central lines
    if origin == True: # origin
        pl.plot([0,0],[np.abs(Y).max()*2,-np.abs(Y).max()*2],'--',color=[0.5,0.5,0.5])
        pl.plot([np.abs(X).max()*2,-np.abs(X).max()*2],[0,0],'--',color=[0.5,0.5,0.5])
    # points
    pl.plot(X,Y,format,**args)
    # names
    
    if names is not None:
        font = {'family' : 'normal',
                'size'   : 12,
                'color'  : (0.7,0.4,0.4)}
        
        
        his = np.histogram2d(X, Y, bins=10)

        
        
        r = []
        for x,y in zip(X,Y):
            r.append((x**2+y**2)**0.5)
            
        #scale = 1-0.8

        max_h = np.percentile(his[0].flatten(), scale*100)
        
        #print his[0].flatten().sum()
        #print his[0].flatten()[np.argwhere(his[0].flatten()>0)]
        #print np.percentile(his[0].flatten()[np.argwhere(his[0].flatten()>0)],scale*100)
        
        max_h = float(np.percentile(his[0].flatten()[np.argwhere(his[0].flatten()>0)],scale*100))
        max_h = float(np.percentile(his[0].flatten()[np.argwhere(his[0].flatten()>0)],scale*100))
        
        #max_r = max(np.array(r).mean(),np.median(np.array(r)))*(1+1-(scale)
        #max_r = max(np.array(r).mean(),np.median(np.array(r)))
        
        max_r = np.percentile(np.array(r), 100-scale*100)
        
        max_valx = max(X)*scale
        max_valy = max(Y)*scale
        min_valx = min(X)*scale
        min_valy = min(Y)*scale

        for x,y,s in zip(X,Y,names):
            r = (x**2+y**2)**0.5
            
            xx = np.argwhere(his[1]>=x)[0]-1
            yy = np.argwhere(his[2]>=y)[0]-1
            #print xx,yy,
            h = float(his[0][xx,yy])
            #print h,max_h
            #if (extreme==False) or ((x>max_valx or x<min_valx) or (y>max_valy or y<min_valy) or (r>max_r)):
            #if (extreme==False) or ((r>max_r)):
            #print r,max_r,s
            #print h,max_h
            if (extreme==False) or ((h<=max_h) or (r>max_r)):
            #if (extreme==False) or ((r>max_r)):
            
                pl.text(x, y, str(s),fontdict=font)
    
    pl.xlabel('PC%d' % PC1)
    pl.ylabel('PC%d' % PC2)
    
    set_axis(X, Y)

