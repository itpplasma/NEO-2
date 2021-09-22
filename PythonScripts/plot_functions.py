#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 14:52:59 2020

@author: wakatobi
"""
__version__='0.0.5'
import os
import matplotlib.pyplot as plt
import h5py
import numpy as np
from matplotlib.figure import Figure
from scipy.interpolate import interp2d

class Neo2Files:
    """Class of File for Neo2 Runs""" # Abstract class usefull?


    def __init__(self,path):

        self._pathInitial=path
        self.exact_file_names=[]
        self.file_paths=[] ### dict with exact_file_name or list?
        if os.path.isfile(path):
            pass


    def get_files_with_exact_names(self):
        """Looking for Filenames from exact_file_names"""

        for root, dirs, files in os.walk(self._pathInitial):
            for exact_file_name in self.exact_file_names:
                if exact_file_name in files:
                    self.file_paths.append(os.path.join(root,exact_file_name))
                    #or if dict: self.dict(exact_name)=[].append(root)
                    #or "/".join(root,exact_name)

class BootstrapFiles(Neo2Files):

        def __init__(self,path):
            super().__init__(path)
            self.exact_file_names.append("g_vs_lambda.h5")

class output_files():
    """Abstract class of outputfiles for NEO2"""

    """One File per Instance"""### One File per Class or other concept??
    def __init__(self):

        self.filedescription="" # Blablabla
        self.filename="" #final.h5
        self.type="output"
        self.fileformat="" # hdf5, fortan namelist


class boozerfile_hdf5(output_files):
    """File of .bc.h5."""

    def __init__():
        super().__init__()


def plot_g(path,title=''):
    """Plot odd part of distribution function

    Needs Reconstruct Mode from Neo2
    """

    if not os.path.isfile(path):
        return 1
    g_vs_h5=h5py.File(path)
    plt.figure()

    for i in g_vs_h5:
        if i == 'version':
            continue

        lam=g_vs_h5[i]['lambda'].value

        g=g_vs_h5[i]['g'].value

        pos_lambda=lam[lam>0]

        pos_g=g[lam>0]

        neg_lambda=lam[lam<0]

        neg_g=g[lam<0]

        odd_g=(pos_g -neg_g[::-1])/2

        odd_lambda=(pos_lambda - neg_lambda[::-1])/2


        odd_eta=1-np.power(pos_lambda,2)


        plt.plot(odd_eta,odd_g)
    plt.xlim(0.8,1)
    plt.xlabel(r'$\eta$')
    plt.ylabel(r'$g_1^a$')
    plt.title(title)

    g_vs_h5.close()

class Plotfunctions():

    def __init__(self,path2bch5):

        if os.path.isfile(path2bch5):
            self.path2bch5=path2bch5
        else:
            raise IOError('Must be a file')
        self.figures=[]

    def add_fig(self):
        self.figures.append(plt.figure())

    def plot_boozer(self,new_figure=True):
        """Plots different representations of boozerfile

        according to the choosen surface
        """

        bch5=h5py.File(self.path2bch5,'r')
        theta_n=bch5['theta_n'].value
        theta=np.array(bch5['theta'])
        phi_n=bch5['phi_n'].value
        phi=np.array(bch5['phi'])

        bmod=np.array(bch5['bmod'])
        bmod=bmod.reshape(phi_n+1,theta_n+1)
        bch5.close()

        if new_figure:
            self.add_fig()
        axs=plt.gca()
        magfield=axs.pcolor(phi,theta,bmod.T)
        cbar=plt.colorbar(magfield)
        axs.set_ylabel(r'$\theta$')
        axs.set_xlabel(r'$\varphi$')
        cbar.set_label('B')

    #Alternativ:
    #plt.contour(phi,theta,bmod.T,20)
    #oder
    #plt.contourf(phi,theta,bmod.T,20)

    def plot_boozer_fieldline(self,len_phi=12,phistart=0,thetastart=0,new_figure=True):
        """Plots the magnetic fieldline in inside the surface"""

        if new_figure:
            fig,ax =plt.subplots()
        self.plot_boozer(new_figure=False)


        bch5=h5py.File(self.path2bch5,'r')

        theta=np.array(bch5['theta'])
        phi=np.array(bch5['phi'])

        maxphi=max(phi)
        maxtheta=max(theta)
        iota=bch5['iota'].value
        bch5.close()



        plt.plot(phistart,thetastart,'bo')

############## Fieldline plotting ############################################

        phiend=(len_phi+phistart)
        while phiend>maxphi:
            thetaend=(thetastart+(maxphi-phistart)*iota)# Total Korrekt

            while thetaend>=maxtheta: # Korrekt
                phiend_i=(maxtheta-thetastart)/iota+phistart
                plt.plot([phistart,phiend_i],[thetastart,maxtheta],'r')
                phistart=phiend_i
                thetaend=thetaend-maxtheta
                thetastart=0

            plt.plot([phistart,maxphi],[thetastart,thetaend],'r')
            thetastart=thetaend
            phiend=phiend-maxphi
            phistart=0


        thetaend=thetastart+phiend*iota


        ###### Same as above ###################
        while thetaend>=maxtheta:
            phiend_i=(maxtheta-thetastart)/iota+phistart
            plt.plot([phistart,phiend_i],[thetastart,maxtheta],'r')
            phistart=phiend_i
            thetaend=thetaend-maxtheta
            thetastart=0
        ########################################


        plt.plot([phistart,phiend],[thetastart,thetaend],'r')
        #plt.show()
        #return ax
##############################################################################



    def plot_boozer_fieldline_poi(self,poi,len_phi=12,phistart=0,thetastart=0,new_figure=True):
        """Plots Point of Intererst into the fieldline"""

        plt.figure()
        plt.subplot(211)
        self.plot_boozer_fieldline(len_phi,phistart,thetastart,new_figure=False)



        bch5=h5py.File(self.path2bch5,'r')


        theta_n=bch5['theta_n'].value
        theta=np.array(bch5['theta'])
        phi_n=bch5['phi_n'].value
        phi=np.array(bch5['phi'])
        bmod=np.array(bch5['bmod'])
        bmod=bmod.reshape(phi_n+1,theta_n+1)

        fit=interp2d(phi,theta,bmod.T)
        maxphi=max(phi)
        maxtheta=max(theta)
        iota=bch5['iota'].value
        bch5.close()


        plot_boozer_poi(poi,iota,maxphi,maxtheta,phistart,thetastart)

        plt.subplot(212)

        def get_B_phiS(bs_l,iota,maxphi,maxtheta,phistart,thetastart):

            bx_l=bs_l/np.sqrt(1+iota*iota)
            by_l=bx_l*iota
            bx_l=bx_l+phistart
            by_l=by_l+thetastart
            bz_l=[]
            for i,j in zip(bx_l%maxphi,by_l%maxtheta):
                bz_l.append(fit(i,j))
            bz_l=np.array(bz_l).flatten()
            return(bz_l)

        bs_l=np.linspace(0,len_phi,1000)*np.sqrt(1+iota*iota)

        bz_l=get_B_phiS(bs_l,iota,maxphi,maxtheta,phistart,thetastart)


        plt.plot(bs_l,bz_l,'r')
        bz_p=get_B_phiS(poi,iota,maxphi,maxtheta,phistart,thetastart)
        plt.plot(poi,bz_p,'+')
        for i,j in zip(poi,bz_p):
            plt.annotate('%s'%i,xy=(i,j))


        plt.ylabel(r'$B$')
        plt.xlabel(r'$\varphi_s$')

def plot_boozer_poi(poi,iota,maxphi,maxtheta,phistart=0,thetastart=0):
    """Plots Point of Interest into the boozer representation"""

    poi_phi,poi_theta=get_phi_theta(poi,iota,maxphi,maxtheta)
    poi_phi=poi_phi+phistart
    poi_theta=poi_theta+thetastart

    for i,j,k in zip(poi,poi_phi%maxphi,poi_theta%maxtheta):
        plt.plot(j,k,c='r', marker='o')
    for i,j,k in zip(poi,poi_phi%maxphi,poi_theta%maxtheta):
        plt.annotate('%s'%i,xy=(j,k))
    plt.ylabel(r'$\theta$')
    plt.xlabel(r'$\varphi$')





def get_phi_theta(bs_p,iota,maxphi,maxtheta):
    """Get phi and thetha out of Fieldline and iota"""

    bx_p=bs_p/np.sqrt(1+iota*iota)
    by_p=bx_p*iota
    return(bx_p,by_p)
