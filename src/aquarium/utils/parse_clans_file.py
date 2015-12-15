#!/usr/bin/env python

import re
import numpy as np
from itertools import izip
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
            

def pairwise(iterable):
    # stackoverflow.com/questions/5389507/iterating-over-every-two-elements-in-a-list#5389547
    "s -> (s0,s1), (s2,s3), (s4, s5), ..."
    a = iter(iterable)
    return izip(a, a)

def sixwise(iterable):
    a = iter(iterable)
    return izip(a, a, a, a, a, a)

def joindicts(dicts):
    for d in dicts:
        d_out = {}
        for dd in d:
            d_out.update(dd)
        yield d_out

def is_iterable(l):
    try:
        _ = (e for e in l)
        return True
    except TypeError:
        pass
    return False


def listify(gen):
    # http://argandgahandapandpa.wordpress.com/2009/03/29/python-generator-to-list-decorator/
    # improved by tljm: for non iterable objects it returns list of one element: the object
    def patched(*args, **kwargs):
        obj = gen(*args, **kwargs)
        if is_iterable(obj):
            return list(obj)
        return [obj]

    return patched

def showit(gen):
    def patched(*args, **kwargs):
        obj = gen(*args, **kwargs)
        plt.show()
        return obj
    return patched


class ClansReader(object):

    
    def __init__(self,filename):
        with open(filename) as f:
            lines = f.readlines()
            self.lines = [line.strip() for line in lines]
        self.default_color = [0.,0.,0.]

    @property
    @listify
    def color_by_group(self):
        seqg = list(self.seqgroups)
        for seq in range(self.sequences):
            in_group = False
            for sg in seqg:
                if seq in sg['numbers']:
                    in_group = True
                    yield ([e/255. for e in sg['color']])[:3]
            if not in_group:
                yield self.default_color
    
    
    @property
    def sequences(self):
        sequences_re = re.compile('^sequences=')
        for line in self.lines:
            match = sequences_re.match(line)
            if match:
                return int(line[match.end():])


    @property
    def seq(self):
        seq_re_beg = re.compile('<seq>')
        seq_re_end = re.compile('</seq>')
        for nr,line in enumerate(self.lines):
            if seq_re_beg.match(line):
                seq_out = []
                for line in self.lines[nr+1:]:
                    if seq_re_end.match(line):
                        return [(seqid[1:],seqseq) for seqid,seqseq in pairwise(seq_out)]
                    seq_out.append(line)
                
    @property
    def seqgroups(self):
        seqgroups_re_beg = re.compile('<seqgroups>')
        seqgroups_re_end = re.compile('</seqgroups>')
        for nr,line in enumerate(self.lines):
            if seqgroups_re_beg.match(line):
                seqgroups_out = []
                for line in self.lines[nr+1:]:
                    if seqgroups_re_end.match(line):
                        return joindicts(sixwise(seqgroups_out))
                    key = line.split('=')[0]
                    value = line.split('=')[1]
                    if key not in ['name']:
                        value = [int(v) for v in value.split(';') if len(v) > 0 ]
                    seqgroups_out.append({key:value})

    @property
    def pos(self):
        pos_re_beg = re.compile('<pos>')
        pos_re_end = re.compile('</pos>')
        for nr,line in enumerate(self.lines):
            if pos_re_beg.match(line):
                pos_out = []
                for line in self.lines[nr+1:]:
                    if pos_re_end.match(line):
                        return np.array(pos_out)
                    pos_out.append(map(float,line.split()[1:]))


    @property
    def hsp_matrix(self):
        hsp_re_beg = re.compile('<hsp>')
        hsp_re_end = re.compile('</hsp>')
        nr_of_sequences = self.sequences
        for nr,line in enumerate(self.lines):
            if hsp_re_beg.match(line):
                hsp_out = np.zeros((nr_of_sequences,nr_of_sequences))
                for line in self.lines[nr+1:]:
                    if hsp_re_end.match(line):
                        return hsp_out
                    xy = map(int,line.split(':')[0].split())
                    xy.sort()
                    x,y = xy
                    hsp_out[x,y] = float(line.split(':')[1])

    @property
    def hsp(self):
        hsp_re_beg = re.compile('<hsp>')
        hsp_re_end = re.compile('</hsp>')
        nr_of_sequences = self.sequences
        for nr,line in enumerate(self.lines):
            if hsp_re_beg.match(line):
                hsp_out = []
                for line in self.lines[nr+1:]:
                    if hsp_re_end.match(line):
                        return np.array(hsp_out)
                    xy = map(int,line.split(':')[0].split())
                    xy.sort()
                    x,y = xy
                    hsp_out.append([x,y,float(line.split(':')[1])])


class ClansReaderPlot(ClansReader):

    def __init__(self,*args,**kwargs):
    
        ClansReader.__init__(self,*args,**kwargs)
        
        self.ax = None
        self.fig = None
    
    @showit
    def init_ax(self):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection='3d') 
        
        self.ax.set_axis_bgcolor('none')
        self.ax.axis('off')
        
        self.fig.subplots_adjust(left=0,bottom=0,right=1,top=1)
        self.fig.set_facecolor('w')

        
    @showit
    def scatter(self,colors=None,s=90,**kwargs):
        if colors is None:
            colors = self.color_by_group
        scat = self.ax.scatter3D(self.pos[:,0],
                                 self.pos[:,1],
                                 self.pos[:,2],
                                 c=colors,s=s,zorder=1,
                                 **kwargs)
        scat.set_clip_on(False)

    @showit
    def text(self,labels=None,shift=None,alpha=0.8,colors='k',**kwargs):
        if shift is None:
            shift = (self.pos.max()-self.pos.min())/100
        if labels is None:
            labels = [s[0] for s in self.seq]
        if not isinstance(colors,list):
            colors = [colors]*self.sequences
        for xyz,label,color in zip(self.pos,labels,colors):
            xyz += shift
            self.ax.text3D(xyz[0],
                           xyz[1],
                           xyz[2],
                           label,color=color,
                           **kwargs)

    @showit
    def forces(self,**kwargs):
        hsp = self.hsp
        force = hsp[:,2]

        empty_force = force == 0
        use_force = force > 0
        force = force[use_force]
        
        force = -np.log10(force)
        force -= force.min()
        force /= force.max()
        
        force = force**(1./4)
        
        force_ = np.zeros(hsp[:,2].shape)
        force_[use_force] = force
        #force_ = 1 - force_
        #force_ = hsp[:,2]
        
        pos = self.pos
        for nr,ab in enumerate(hsp[:,:2]):
            a,b = ab
            # empty?
            #if empty_force[nr]:
            #    continue
            f = force_[nr]
            #if f > threshold:
            #    continue
            xyz = pos[[a,b],:]
            self.ax.plot3D(xyz[:,0],
                           xyz[:,1],
                           xyz[:,2],
                           c=[f,f,f],alpha=(1-f)/2.,zorder=-1,
                           **kwargs)

