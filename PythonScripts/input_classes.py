#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  9 17:21:45 2018

@author: wakatobi
"""
### Template from Ipthon Notebook doc_18_07_05

import sys
sys.path.append('/afs/itp.tugraz.at/user/wakatobi/Modul_neo/f90nml-0.21/')
import f90nml
import os
from os.path import *


class Multitest():
    ### Methods should be used but not class in total.
    
    def __init__(self,runpath,runsource='/afs/itp.tugraz.at/user/wakatobi/NEO_2_V0.10/NEO-2-QL/Build5/neo_2.x'):
        self.runpath=runpath
        #self.codesource=None
        #self.runsource=None # Executable
        if isfile(runsource):
            self.runsource=runsource
        
        self.reqfiles={
        'neo2in': 'neo2.in',
        'neoin': 'neo.in'
        }
        self.sourcepaths=dict()
        self._SetSources()
        self._SetSourcePaths() # including Fill required Files from neo Files, which reads neo2.in
        
        
        ##########Info :self.sources=dict()
       
        #self.SourcePaths=dict()
        #self.SetSourcepaths()

    
    
    #### INIT DONE FOR THE FIRST!!!!####
    
    
    
    
    
    def Run_Precomp(self,overwrite=False):
        
        #load_parameter from an old run! 
        self.CheckReqFiles()
        self.CreateLinkFiles(overwrite)
        self.RunInitRuns()
        #self.gen_condorinput()
        #self.run_condor()
        
        
        
    def RunInitRuns(self):
        ###Use Routine from andfmar
        
        if self._neo_nml['multi_spec']['isw_multispecies_init'] != 1:
            raise InputError('isw_multispecies_init is not 1')
        else:
            curdir=os.getcwd()
            print('curdir = ',curdir)
            print('runpath = ',self.runpath)
            os.chdir(self.runpath)
            print('Starting initrun')
            print(os.popen("./neo_2.x").read())
            os.chdir(curdir)
        

    def CreateLinkFiles(self,overwrite=False,link=True):
        if self.runpath==None:
            raise IOError('runpath is missing')
        os.makedirs(self.runpath,exist_ok=True)
        
        for i in self.sourcepaths.values():
            print(i)
            destfile=join(self.runpath,basename(i))
                
            if isfile(destfile):
                if overwrite:
                    print(destfile, ' is a file and is replaced')
                    os.remove(destfile)
                else:
                    print(destfile, ' is not updated') ### Maybe check if it is the same file!!
                    continue
            
            os.symlink(i,destfile)
        
        
        
        i=self.runsource
        print('runsource = ', i)
        destfile=join(self.runpath,basename(i))
        print('destfile = ', destfile)
        if isfile(destfile):
            if overwrite:
                print(destfile, ' is a file and is replaced')
                os.remove(destfile)
            else:
                print(destfile, ' is not updated') ### Maybe check if it is the same file!!
                return
            
        os.symlink(i,destfile)             

        
        
        
    def CheckReqFiles(self):
        
        if set(self.reqfiles) != set(self.sourcepaths):
            raise ValueError('Not all paths of required Files defined')
        
        for i,j in self.sourcepaths.items():
            if not exists(j):
                raise RuntimeError(j,' is not a File')
                
            
        else:
            print('All required Files are existing')
                

    
        
    def SetNeo2file(self,neo2file): ## Method to set the path of the neo2 file and make all check(Setter of neo2.in)
        if os.path.isfile(neo2file):
            self.sourcepaths['neo2in']=os.path.abspath(neo2file)
            
        if isdir(neo2file):
            self.sourcepath['neoin']=join(abspath(neo2file),'neo.in')
            self.sourcepaths['neo2in']=join(abspath(neo2file),'neo2.in')   
            
            
            
            
            
        
        
        
####### Implementend Functions: #########       
        
        
    def _FillReqFiles(self): # Get required File paths from neo(2).in Files
        
        self._read_neo2in()
        try:
            self.reqfiles['multispec']=self._neo_nml['multi_spec']['fname_multispec_in']
        except KeyError:
            print('The neo2.in File has no multispec parameter')
            return 
        
        
        self.reqfiles['in_file_pert']=self._neo_nml['ntv_input']['in_file_pert']
        
       
        
        try:
            neo=self.sourcepaths['neoin']
            ## There must be a better catch if neoin is not a File
        except:
            print('neo.in File is missing')
            return
        
        with open(neo,'r') as f:
            for line in f:
                arg = line.split()
                if 'in_file' in arg:
                    self.reqfiles['in_file_axi']=arg[0]
                    #print(arg[0])
                    break 
    # Ensure that it is running in the correct mode
    
    def _SetSources(self): # name has to be better chosen:
        
        ##### Uniqifi path if only a oberpath is given
        self.sources=dict([('codesource','/proj/plasma/Neo2/MODULAR/'),
                                ('pert_path','/temp/andfmar/Boozer_files_perturbation_field/ASDEX_U/32169/'),                            
                                ('singlerunsource','/temp/wakatobi/Startordner_Neo2/')])
        
        
    def _SetSourcePaths(self,overwrite=False): # name has to be better choosen
        
        try:
            files=os.listdir(self.sources['singlerunsource'])
        except:
            print('sources not correct set')
            return
        
        for fdisc,fnames in self.reqfiles.items():
            
            if fnames in files:
                if fnames in self.sourcepaths and overwrite==False:
                    #print(fnames, ' is already set to path: ' self.sourcepaths[fnames])
                    continue
                self.sourcepaths[fdisc]=join(self.sources['singlerunsource'],fnames)
                
        self._FillReqFiles()
        #read_Neo2File
        
        
        ### DOTO Implement method to iterate over all sources, otherwise problems are occuring
        if set(self.reqfiles) == set(self.sourcepaths):
            return
        
        try:
            files=os.listdir(self.sources['pert_path'])
        except:
            print('sources2 not correct set')
            return
        
        for fdisc,fnames in self.reqfiles.items():
            
            if fnames in files:
                if fnames in self.sourcepaths:
                    #print(fnames, ' is already set to path: ' self.sourcepaths[fnames])
                    continue
                self.sourcepaths[fdisc]=join(self.sources['pert_path'],fnames)        

    
    def _read_neo2in(self):
        try:
            self._neo_nml=f90nml.read(self.sourcepaths['neo2in'])
        except:
            print('Couldn\'t read neo2.in')
            return
