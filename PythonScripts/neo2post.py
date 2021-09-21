#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 16:43:38 2020

@author: wakatobi
"""

import matplotlib.pyplot as plt
import h5py



def myplot(*args,ax=None,def_x=None,scalex=True, scaley=True, **kwargs):
    ''' Function for plotting  values, possible with a default x_axis

    args and def_x should be Numpy Arrays
    '''

    if ax==None:
        fig=plt.figure()
        ax=fig.add_subplot(111)
        ax.set_xlabel('boozer_s')

    if len(args)==1:
        #ax.set_xlabel('boozer_s')
        if len(args[0].shape)>2:
            pass
            #return ax.plot()
        else:
            return ax.plot(def_x,*args,scalex=scalex,scaley=scaley,**kwargs)

    else:
        raise NotImplementedError

plot_h5py_1d=myplot
plot_h5py_stack_1d=myplot
def plot_h5py_Dataset(dataset,*args,**kwargs):
    if dataset.ndim == 1:
        plot_h5py_1d(dataset,*args,**kwargs)
    elif dataset.ndim == 2:
        if dataset.shape[1]<9:
            plot_h5py_stack_1d(dataset,*args,**kwargs)
        else:
            pass
            #plot_h5py_2d(dataset,*args,**kwargs)
    else:
        pass
        #raise NotImplementedError("to many dimensions for a Line Plot")


class Neo2Plot():
    '''Plotting all types of NEO2Files

    '''
    def_x='boozer_s'
    def __init__(self,file):
        self.NA_list=['D11_NA_Dpl', 'D12_NA_Dpl', 'D21_NA_Dpl', 'D22_NA_Dpl']
        self.file=file ## So far file must be h5py


    def plot(self,*args,**kwargs):
        fig=plt.figure()
        for ind, i in enumerate(args):
                ax=fig.add_subplot(4,1,ind+1)
                plot_h5py_Dataset(self.file[i].value,ax=ax,def_x=self.file[self.def_x].value,**kwargs)
                ax.grid(True)
                ax.set_ylabel(i)
        fig.tight_layout()

    def plot_NA(self):
        fig=plt.figure()

        for ind,i in enumerate(['D11_NA_Dpl', 'D12_NA_Dpl', 'D21_NA_Dpl', 'D22_NA_Dpl']):
                ax=fig.add_subplot(4,1,ind+1)
                plot_h5py_Dataset(self.file[i].value,ax=ax,def_x=self.file[self.def_x].value)
                ax.grid(True)
                ax.set_ylabel(i)
        fig.tight_layout()
    def __dir__(self):
        '''Function for displaying Parameter as attribute'''
        return(dir(Neo2Plot)+list(self.final))

    def __getattr__(self,name):
        '''Function for plotting one Parameter as attribute'''
        if name.startswith('_'):
            raise AttributeError(name)
        if name not in list(self.file):
            raise AttributeError(name)
        self.plot(name)
    def __repr__(self):
        return "\n".join(list(self.file))
    def _ipython_key_completions_(self):
        '''Auto complete possible plot parameter'''
        return list(self.file)
    def __getitem__(self,key):
        '''Also subscribing desired parameter to plot  is possible'''
        return self.plot(key)

class Multirun():
    '''Generel class for Neo2 run collections'''
    def __init__(self,initdir):

        self.initdir=initdir
        if h5py.is_hdf5(initdir):
            self.scanfile=initdir
        self.singleruns=[]
        self.postprocessed=[]
        self.otherfiles=[]
        self.otherdirs=[]



class RadialScan(Multirun):
    '''Multiple Neo2runs with different radial component'''


    def_x='boozer_s'
    def __init__(self,initdir):
        super().__init__(initdir)
        self.NA_list=['D11_NA_Dpl', 'D12_NA_Dpl', 'D21_NA_Dpl', 'D22_NA_Dpl']

        self.finalhdf5=h5py.File(self.scanfile)
        self.fig=Neo2Plot(self.finalhdf5)
