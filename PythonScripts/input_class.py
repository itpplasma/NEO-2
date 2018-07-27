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
import yaml
from IPython import display

class Multitest():
    ### Methods should be used but not class in total.
    
    def __init__(self,runpath):
        self.runpath=runpath
        #self.codesource=None
        #self.runsource=None # Executable
    
        self._reqfiles={
        'neo2in': 'neo2.in',
        'neoin': 'neo.in'
        }
        self._Load_default()
        self._sourcepaths=dict()
        #self._SetSources()
        self._SetSourcePaths() # including Fill required Files from neo Files, which reads neo2.in
        
        
        ##########Info :self._sources=dict()
       
        #self.SourcePaths=dict()
        #self.SetSourcepaths()

    
    
    #### INIT DONE FOR THE FIRST!!!!####
    
    
    def _Load_default(self):
        try:
            self.__sourcedir=__file__
        except NameError:
            self.__sourcedir=''
        
        self.__realpath=os.path.realpath(self.__sourcedir)
        if os.path.isfile(self.__realpath):
            path = os.path.dirname(os.path.realpath(self.__sourcedir))
        elif os.path.isdir(self.__realpath):
            path=self.__realpath
    
        yamlpath = os.path.join(path, 'input_default.yaml')
        print(path)
        try:
          with open(yamlpath) as f:
              data=yaml.load(f)
        except:
            print('Couldn\'t read default values')
            return
        
        self._runsource=data.get('runsource')
        self._reqfiles=data.get('reqfiles')
        self._sources=data.get('sources')
    
    
    
    def Run_Precomp(self,overwrite=False):
        
        #load_parameter from an old run! 
        self._CheckReqFiles()
        self.CreateFiles(overwrite)
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
        

    def CreateFiles(self,overwrite=False,link=True):
        if self.runpath==None:
            raise IOError('runpath is missing')
        os.makedirs(self.runpath,exist_ok=True)
        
        for i in self._sourcepaths.values():
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
        
        
        
        i=self._runsource
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

        
        
        
    def _CheckReqFiles(self):
        
        if set(self._reqfiles) != set(self._sourcepaths):
            raise ValueError('Not all paths of required Files defined')
        
        for i,j in self._sourcepaths.items():
            if not exists(j):
                raise RuntimeError(j,' is not a File')
                
            
        else:
            print('All required Files are existing')
                

    
        
    def SetNeo2file(self,neo2file): ## Method to set the path of the neo2 file and make all check(Setter of neo2.in)
        if os.path.isfile(neo2file):
            self._sourcepaths['neo2in']=os.path.abspath(neo2file)
            
        if isdir(neo2file):
            self.sourcepath['neoin']=join(abspath(neo2file),'neo.in')
            self._sourcepaths['neo2in']=join(abspath(neo2file),'neo2.in')   
            
            
            
            
            
        
        
        
####### Implementend Functions: #########       
        
        
    def _FillReqFiles(self): # Get required File paths from neo(2).in Files
        
        self._read_neo2in()
        try:
            self._reqfiles['multispec']=self._neo_nml['multi_spec']['fname_multispec_in']
        except KeyError:
            print('The neo2.in File has no multispec parameter')
            return 
        
        
        self._reqfiles['in_file_pert']=self._neo_nml['ntv_input']['in_file_pert']
        
       
        
        try:
            neo=self._sourcepaths['neoin']
            ## There must be a better catch if neoin is not a File
        except:
            print('neo.in File is missing')
            return
        
        with open(neo,'r') as f:
            for line in f:
                arg = line.split()
                if 'in_file' in arg:
                    self._reqfiles['in_file_axi']=arg[0]
                    #print(arg[0])
                    break 
    # Ensure that it is running in the correct mode
    
    #def _SetSources(self): # name has to be better chosen:
        
        ##### Uniqifi path if only a oberpath is given
       # self._sources=dict([('codesource','/proj/plasma/Neo2/MODULAR/'),
                           #     ('pert_path','/temp/andfmar/Boozer_files_perturbation_field/ASDEX_U/32169/'),                            
                #                ('singlerunsource','/temp/wakatobi/Startordner_Neo2/')])
        
        
    def _SetSourcePaths(self,overwrite=False): # name has to be better choosen
        
        try:
            files=os.listdir(self._sources['singlerunsource'])
        except:
            print('sources not correct set')
            return
        
        for fdisc,fnames in self._reqfiles.items():
            
            if fnames in files:
                if fnames in self._sourcepaths and overwrite==False:
                    #print(fnames, ' is already set to path: ' self._sourcepaths[fnames])
                    continue
                self._sourcepaths[fdisc]=join(self._sources['singlerunsource'],fnames)
                
        self._FillReqFiles()
        #read_Neo2File
        
        
        ### DOTO Implement method to iterate over all sources, otherwise problems are occuring
        if set(self._reqfiles) == set(self._sourcepaths):
            return
        
        try:
            files=os.listdir(self._sources['pert_path'])
        except:
            print('sources2 not correct set')
            return
        
        for fdisc,fnames in self._reqfiles.items():
            
            if fnames in files:
                if fnames in self._sourcepaths:
                    #print(fnames, ' is already set to path: ' self._sourcepaths[fnames])
                    continue
                self._sourcepaths[fdisc]=join(self._sources['pert_path'],fnames)        

    
    def _read_neo2in(self):
        try:
            self._neo_nml=f90nml.read(self._sourcepaths['neo2in'])
        except:
            print('Couldn\'t read neo2.in')
            return






class Neo2File(object):
    
    def __init__(self,neo2path):
        self._neo2nml=f90nml.read(neo2path)
        self._neo2dict=dict()
        self._nmltodict()
    
    
    def _nmltodict(self):
        for i in self:
            self._neo2dict[i]=self[i]
    
    def chval(self,**kwargs):
        
        print(self['plotting'])
        for par in kwargs:
            
            match=False
            for i in self._neo2nml: # Iterating through namelists
                if par in self._neo2nml[i]:  # if parameter in one namelist 
                    if self._neo2nml[i][par]==kwargs[par]:
                        print(k, ' in namelist', i, ' was already : ', kwargs[par])
                    else:
                        print('before change: ',par,'=',self._neo2nml[i][par])
                        self._neo2nml[i][par]=kwargs[par]
                        print(par, ' in namelist', i, ' changed to : ', kwargs[par])
                    match=True
                    break
                if match:
                    break
            else:
                print(par, ' is not in the neo2.in File')
    
    
    
    def write(self):
        return self._neo2nml.write()
    
    
    
    def __repr__(self):
        display.display_pretty(self._neo2dict)
    
    
    #def __repr__(self):
     #   lines = ['{']
      #  for key, value in self._neo2dict.items():
       #     lines.append('{}:{}'.format(key, value))
        #lines.append(['}'])
        
       # return '\n'.join(lines)
    
    
    #def __getattr__(self, attr):
      #  if attr=='write':
       #     return getattr(self._neo2nml,attr)

    #def __dir__(self):
     #   print('__dir__ is called')
    #   return ['write']
        #return self._neo2nml.__dir__()
       

    
    #def __setattr__(self, attr, val):
       # print('attr = ',attr)
      #  if attr == '_neo2nml':
       #     object.__setattr__(self, attr, val)

       # return setattr(self._neo2nml, attr, val)
    
    def __getitem__(self, key):
        
        ## Unbedingt Kontrollen einführen  
        if key in self._neo2nml:
            print('Please adress your desired parameter directly')
            return None
        else:
            for i in self._neo2nml:
                if key in self._neo2nml[i]:
                    #print('parameter: ',key,'found in namelist:',i )
                    return self._neo2nml[i][key]
        val = self._neo2nml.__getitem__(key)
        if isinstance(val,f90nml.namelist.Namelist):
            valout='namelist of key'
            
        else:
            valout=val
        print('getter: called: ',key,' with value:', valout)
        return val

    def __setitem__(self, par, val):
               
        if par in self._neo2nml:
            print('Please adress your desired parameter directly')
            return None
        else: 
            for i in self._neo2nml:
                if par in self._neo2nml[i]:
                    if not self._checktype(self._neo2nml[i][par],val):
                        print('wrong type of parameter', par,type(val),'-- old type was:', type(self._neo2nml[i][par]))
                        return
                
                    if self._neo2nml[i][par]==val:
                        print(par, ' in namelist', i, ' was already : ', val)
                        return
                    else:
                        print('before change: ',par,'=',self._neo2nml[i][par])
                        self._neo2nml[i][par]=val
                        print(par, ' in namelist', i, ' changed to : ', val)
                        return
        # Auch unbedingt Kotrollen einführen
        #print('setter:',key,'is set to :', val)
        #self._neo2nml.__setitem__(key, val)
        #write File when set.
        
    def _checktype(self,val1,val2):
        if type(val1)==type(val2):
            return True
        else: 
            return False
        

    def __iter__(self):
       # def chain(*iterables):
            # chain('ABC', 'DEF') --> A B C D E F
        for it in self._neo2nml:
            for element in self._neo2nml[it]:
                yield element
    
        


    

    
    #def update(self, *args, **kwargs):
     #   print('updateoe', args, kwargs)
      #  for k, v in f90nml.namelist.Namelist(*args, **kwargs).iteritems():
       #     self[k] = v
    
   # def __contains__(self, key):
    #    #print('contains called for ',key)
    #    return self._neo2nml.__contains__(key)

















