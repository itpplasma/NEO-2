# -*- coding: utf-8 -*-
"""
Created on Thu May 11 18:18:49 2017

@author: wakatobi
"""
import matplotlib.pyplot as plt
import h5py
import numpy as np
import os
from collections import defaultdict


class Output():

    def __init__(self, path):
        self._pathInitial = path

        self._plot_dict = {}
        self._paradict = {}
        self.list_para = []
        self._h5fsingle = 'neo2_multispecies_out.h5'
        self._h5f = 'final_neo2_multispecies_out.h5'
        self._path = ''
        self._begining()
        #self._SingleRunCheck()

    def _begining(self):    # Check if collected final h5 files is available in
                            # entered path, or in subdirectories.
        self.finalh5_Paths = []
        self.Singleh5_Paths = []
        for root, dirs, files in os.walk(self._pathInitial):
            if self._h5f in files:
                self.finalh5_Paths.append(root)

            if self._h5fsingle in files:
                self.Singleh5_Paths.append(root)
        if len(self.finalh5_Paths) == 0:
            #print('no file with name ' + self._h5f + ' found in ' + self._pathInitial)
            pass
        elif len(self.finalh5_Paths) == 1:
            self._path = self.finalh5_Paths[0]
            self._fill_plot_dict()
            print('one ' + self._h5f + ' File found. Change path to ' + self._pathInitial)
        else:
            print('more than one '+self._h5f + ' File found')

    def _SingleRunCheck(self):
        if isinstance(self._pathInitial, str):
            self._path = self._pathInitial
            self._fill_plot_dict()


    def _collect_plot_data(self):
        pass

    def _check_path(self):

       pass

    def _fill_multirun_dict(self, paths, values, para2check = ''): # values = list of parameter which should be compared
        self._multirunLegend = defaultdict(list)
        self._multirunData = defaultdict(list)
        for i2 in paths:
            print(i2 + ' = path')
            self._ProfileRun = {}
            self._ProfileRun['info'] = []
            self._ProfileRun['x-Achse'] = []
            self._ProfileRun['y-Achse'] = []

            ##
            f = h5py.File(os.path.join(i2, 'final_neo2_multispecies_out.h5'))
            for booz in f:
                if booz == 'version':
                    # print('version')
                    continue
                for para in f[booz]:
                    if len(f[booz][para].shape) >= 2:
                        continue
                    if para not in values:
                        # print(para + ' not in values')
                        continue
                    if para not in self._ProfileRun:
                        self._ProfileRun[para] = []
                        for i in range(f[booz][para].shape[0]):
                            self._multirunLegend[para].append(i2 + ' ['
                                                              + str(i) + ']')

                    self._ProfileRun[para].append(f[booz][para].value)
                    # if para not in self._ProfileRun['info']:
                        # self._ProfileRun['info'].append(para)
                        # self._ProfileRun['x-Achse'].append(para)
                        # self._ProfileRun['x-Achse'][para]=[]
                        # self._ProfileRun['y-Achse'].append(para)
                        # self._ProfileRun['y-Achse'][para]=[]
            # print(self._ProfileRun2)
            for i in self._ProfileRun:
                if i not in values:
                    continue
                self._multirunData[i].append(np.vstack(self._ProfileRun[i]))
                # self._ProfileRun


    def _fill_Single_multirun_dict(self, paths, values):
        self._SingleMultirunLegend = defaultdict(list)
        self._SingleMultirunData = defaultdict(list)
        for i in paths:
            #self._Run = {}

            f = h5py.File(os.path.join(i, 'neo2_multispecies_out.h5'))
            for para in f:

                if para not in values:
                    continue
                if len(f[para].shape) >= 2:
                    print(para+ ' is not a scalar')
                    continue
                self._SingleMultirunData[para].append(f[para].value)
                self._SingleMultirunLegend[para].append(' ['+str(i) + ']')






    def _fill_plot_dict(self):
        f = h5py.File(os.path.join(self._pathInitial,
                                   'final_neo2_multispecies_out.h5'))
        for i in f:
            if i == 'version':
                continue

            for para in f[i]:
                if len(f[i][para].shape) <2:  #Check only 1D and 2D Data

                    if para not in self._plot_dict:
                        self._plot_dict[para]=[]
                    self._plot_dict[para].append(f[i][para].value)
        for i in self._plot_dict:
            self.list_para.append(i)
            self._paradict[i]=np.vstack(self._plot_dict[i])
            if self._paradict[i].shape[-1] > 10:
                self._paradict[i]=np.vstack(self._plot_dict[i]).T

    def plot(self,y_Achse,x_Achse='boozer_s',linstyle='') :
        if x_Achse in self._paradict and y_Achse in self._paradict:
            plt.figure()
            plt.plot(self._paradict[x_Achse],self._paradict[y_Achse],linstyle)
            plt.xlabel(x_Achse)
            plt.ylabel(y_Achse)
        else:
            print(x_Achse,y_Achse,' is not in available. Please check Parameterlist')
    def ambipolarplot(self,linstyle=''):
        ambi=np.sum(self._paradict['Gamma_AX_spec']*self._paradict['z_spec'],1)
        plt.figure()
        plt.plot(self._paradict['boozer_s'],ambi/self._paradict['Gamma_AX_spec'][:,0],linstyle)

        plt.xlabel('boozer_s')
        plt.ylabel('Ambipolar_rel')
    def collect_multirunData(self):
        pass
