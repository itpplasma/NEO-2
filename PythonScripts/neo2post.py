#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 16:43:38 2020

@author: wakatobi
"""

import matplotlib.pyplot as plt
import h5py
import os
import numpy as np


def myplot(*args,ax=None,def_x=None,scalex=True, scaley=True, **kwargs):
    ''' Function for plotting  values, possible with a default x_axis

    args and def_x should be Numpy Arrays
    '''

    if ax:
        plt.sca(ax)
    if len(args)==1:
        #ax.set_xlabel('boozer_s')
        if len(args[0].shape)>2:
            pass
            #return ax.plot()
        else:
            return plt.plot(def_x,*args,scalex=scalex,scaley=scaley,**kwargs)

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
            if 'second_dimension_index' in kwargs:
                second_dimension_index=int(kwargs.pop('second_dimension_index'))
                plot_h5py_stack_1d(dataset[:,second_dimension_index],*args,**kwargs)
            else:
                pass
    else:
        pass
        #raise NotImplementedError("to many dimensions for a Line Plot")


class Neo2Plot():
    '''Plotting all types of NEO2Files

    '''

    def __init__(self,file_,def_x=''):
        self.NA_list=['D11_NA_Dpl', 'D12_NA_Dpl', 'D21_NA_Dpl', 'D22_NA_Dpl']
        if not isinstance(file_,(h5py.File,h5py.Group)):
            raise AttributeError('File must be h5py')
        self.file=file_ ## So far file must be h5py
        if not def_x:
            self.def_x='boozer_s'
        else:
            self.def_x=def_x

        self._valid_keys=list()
        self._get_valid_keys()

    def _get_valid_keys(self):
        x_shape=self.file[self.def_x].shape
        for i in self.file:
            try:
                if self.file[i].shape[0] == x_shape[0]:
                    self._valid_keys.append(i)
            except:
                pass


    def plot(self,*args,**kwargs):

        if 'def_x' in kwargs:
            def_x=kwargs.pop('def_x')
        else:
            def_x=self.def_x
        if 'label' in kwargs:
            label=kwargs.pop('label')
        else:
            label=''

        ## iter_ax  =  0.. No axes
        ##          =  1.. One axes
        ##          =  2.. Multiple axes

        if 'ax' in kwargs:
            ax=kwargs.pop('ax')
            if isinstance(ax,(list)):
                if len(ax)==1:
                    iter_ax=1
                    ax=ax[0]
                elif len(ax)==len(args):
                    iter_ax=2
                else:
                    raise IOError('Number of axes does not match with arguments')
            else:
                iter_ax=1
        else:
            iter_ax=0


        for arg in args:
            arglabel=arg + '/' + label

            if iter_ax==0:
                plot_h5py_Dataset(self.file[arg].value,def_x=self.file[def_x].value,label=arglabel,**kwargs)
            elif iter_ax==1:
                if def_x=='lambda':
                    ax.set_xlabel(r'$\lambda$')
                else:
                    ax.set_xlabel(def_x)
                plot_h5py_Dataset(self.file[arg].value,ax=ax,
                                  def_x=self.file[def_x].value,label=arglabel,**kwargs)
            elif iter_ax==2:
                i_ax=ax.pop()
                if def_x=='lambda':
                    i_ax.set_xlabel(r'$\lambda$')
                else:
                    i_ax.set_xlabel(def_x)
                i_ax.set_ylabel(arg)
                plot_h5py_Dataset(self.file[arg].value,ax=i_ax,
                                  def_x=self.file[def_x].value,label=label,**kwargs)


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
        return(dir(Neo2Plot)+self._valid_keys)

    def __getattr__(self,name):
        '''Function for plotting one Parameter as attribute'''
        if name.startswith('_'):
            raise AttributeError(name)
        if name not in list(self.file):
            raise AttributeError(name)
        self.plot(name)
    #def __repr__(self):
        #return "\n".join(list(self.file))
    def _ipython_key_completions_(self):
        '''Auto complete possible plot parameter'''
        return self._valid_keys
    def __getitem__(self,key):
        '''Also subscribing desired parameter to plot  is possible'''
        return self.plot(key)

    def __call__(self,*args,**kwargs):
        return self.plot(*args,**kwargs)

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

class Bootstrap(Multirun):
    '''Plot up to distribtion function'''

    def __init__(self,initdir):
        files=os.listdir(initdir)
        if 'final.h5' not in files:
            raise IOError('This is not an reconstruction run')
        super().__init__(initdir)

    def interactive_plot(self):
        pass


class RadialScan(Multirun):
    '''Multiple Neo2runs with different radial component'''


    def_x='boozer_s'
    def __init__(self,initdir):
        super().__init__(initdir)
        self.NA_list=['D11_NA_Dpl', 'D12_NA_Dpl', 'D21_NA_Dpl', 'D22_NA_Dpl']

        self.finalhdf5=h5py.File(self.scanfile)
        self.plot=Neo2Plot(self.finalhdf5,def_x=self.def_x)


class MagneticsPlot():
    '''Plot magnetics.h5 after Neo2 Run

    rundir must be an absolute path,
    plotdir can be a relative path to the rundir or an absolute path
    '''

    def __init__(self,rundir,plotdir,templatedir=''):

        if os.path.isdir(rundir):
            self.rundir=rundir
        else:
            raise IOError('rundir must a valid path')

        self.plotdir=os.path.join(rundir,plotdir)
        # if plotdir is absolute, join only outputs plotdir

        self.poi=[]
        self._read_magnetics()

    def _read_magnetics(self):

        magneticsh5=h5py.File(os.path.join(self.rundir,'magnetics.h5'),'r')
        self.s_0,self.phi_0,self.theta_0 = np.array(magneticsh5['fieldline']['1']['xstart'])

        self.iota=magneticsh5['surface']['1']['aiota'].value
        props=[int(i) for i in magneticsh5['fieldpropagator']]
        props.sort()
        bhat_line=[]
        phi_bhat_line=[]
        for prop in props:

            bhat=np.array(magneticsh5['fieldpropagator'][str(prop)]['bhat'])
            bhat_line.append(1/bhat)
            phi=np.array(magneticsh5['fieldpropagator'][str(prop)]['x2'])
            phi_bhat_line.append(phi)

        self.phi_bhat_line_s=np.concatenate(phi_bhat_line[:-1])
        self.bhat_line_s= np.concatenate(bhat_line[:-1])
        self.phi_min=self.phi_bhat_line_s[0]
        self.phi_max=self.phi_bhat_line_s[-1]
        self.phi_begin=self.phi_min
        self.phi_end=self.phi_max
        magneticsh5.close()

    def plot_magnetics(self,begin=None,end=None):

        phi_bhat_line_s=self.phi_bhat_line_s
        bhat_line_s= self.bhat_line_s

        if begin:
            logical_begin=phi_bhat_line_s>begin
            bhat_line_s= bhat_line_s[logical_begin]
            phi_bhat_line_s=phi_bhat_line_s[logical_begin]

        if end:
            logical_end=phi_bhat_line_s<end
            phi_bhat_line_s=phi_bhat_line_s[logical_end]
            bhat_line_s= bhat_line_s[logical_end]

        plt.plot(phi_bhat_line_s,bhat_line_s)
        plt.xlabel(r'$\varphi_s$')
        plt.ylabel(r'$1/\hat{B}$')


    def _plot_singlepoi(self,point,add_point=True):

        point_ind=np.where(self.phi_bhat_line_s-point>0)[0][0]
        plt.plot(self.phi_bhat_line_s[point_ind],self.bhat_line_s[point_ind],marker='o')

        if add_point:
            theta_point = self.theta_0 + self.iota *(self.phi_bhat_line_s[point_ind] - self.phi_0);
            phi_point   = self.phi_bhat_line_s[point_ind]
            theta_point = theta_point%(2*np.pi)
            phi_point = phi_point%(2*np.pi)
            self.poi.append((self.s_0,theta_point,phi_point))
    def plot_poi(self,point,new_points=True):
        if new_points:
            self.poi=[]
        self.plot_magnetics(begin=self.phi_begin,end=self.phi_end)
        if isinstance(point,(list,tuple)):
            for i in point:
                self._plot_singlepoi(i)
        else:
            self._plot_singlepoi(point)

    def write_poi(self,overwrite=False):

        writepath=os.path.join(self.plotdir,'g_vs_lambda.in')
        os.makedirs(self.plotdir,exist_ok=True)
        if overwrite:
            os.remove(writepath)
        mode='x'
        with open(writepath,mode) as f:
            f.write('{:>18s}{:>18s}{:>18s}{:>7s}\n'.format('s', 'theta', 'phi', 'tag'))
            for i,j in enumerate(self.poi):
                f.write('{:18.8e}{:18.8e}{:18.8e}{:>7s}\n'.format(*j,'p'+str(i)))
