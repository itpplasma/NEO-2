#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on ?

@author: wakatobi
"""
class Fluxsurface():

  def __init__(self, rundir='', templatepath='', new_dir=False):

    self._templatepath=templatepath
    self._new_dir=new_dir
    self.filenames={
            'surfaces_in':'create_surfaces.in',
            'spitzerinterface':'spitzerinterface.in',
            'profile':'profiles.dat',
            'neoin':'neo.in',
            'surface_exe':'create_surface.x'} # Create Surf.py, #neo2.in
    self.filepaths=dict(
            neoin='/temp/gernot_k/Neo2/AG/Interface/Productive/2016_02_w7x-m24li/neo.in',
            surface_exe='/proj/plasma/Neo2/Interface/Create_Surfaces/Build/'
            'create_surfaces.x',
            spitzerinterface='/proj/plasma/Neo2/Interface/Create_Surfaces/'
            'create_surfaces.in',
            profile='/proj/plasma/Neo2/Interface/Profiles/w7x-m111-b3-i1/profiles.dat')# Create Surf.py

    if os.path.isdir(rundir) or new_dir :
      self._rundir=rundir
    else:
      raise IOError('rundir must a valid path')

    self.get_filenames()
    self.get_filepaths(path='/temp/gernot_k/Neo2/AG/Interface/Productive/2016_02_w7x-m24li/',
                               replace=False) ## Fühlt nur missende auf,(in_file_axi)
    self._createfiles()
    if new_dir:
      self._run_surface()

  def get_filenames(self):
    # Check of filename inside surfaces.in has to be done.

    try:
      neo=self.filepaths['neoin']
      ## There must be a better catch if neoin is not a File
    except:
      print('neo.in File is missing')
      return

    with open(neo,'r') as f:
      for line in f:
        arg = line.split()
        if 'in_file' in arg:
          self.filenames['in_file_axi']=arg[0]
          #print(arg[0])
          break


  def get_filepaths(self, path='', replace=True):

    #TODO reset option if path is changing
    if not path:
      path=self._templatepath

    try:
      files=os.listdir(path)
    except:
      print('path is not set correctly')
      return

    for file,filename in self.filenames.items():

      if filename in files:
        if file in self.filepaths and replace==False:
          print(file, ' is already set to path: ', self.filepaths[file])
          continue
    ##TODO find orginal path to files not only links
        self.filepaths[file]=os.path.join(path,filename)

    ### DOTO Implement method to iterate over all sources, otherwise problems are occuring
    if set(self.filenames) == set(self.filepaths):
      print('filenames and filepaths agree')
    else:
      print('filenames and filepaths do not agree')


  def _run_surface(self,show=True):
    curdir=os.getcwd()
    os.chdir(self._rundir)

    os.chdir(curdir)


  def _createfiles(self, rundir='', overwrite=False, link=False):

    if not rundir:
        rundir=self._rundir
    if self._new_dir:
        os.mkdir(self._rundir)
        self._new_dir=False
    #self._checkreqfiles() #ToDO Make own _checkreqfiles()

    #os.makedirs(newplotdir,exist_ok=True)

    for i,j in self.filepaths.items():

      destfile=os.path.join(rundir,os.path.basename(j))

      if os.path.isfile(destfile):
        if overwrite:
          print(destfile, ' is a file and is replaced')
          os.remove(destfile)
        else:
          print(destfile, ' is not updated') ### TODO Maybe check if it is the same file!!
          continue
      if link:
        os.symlink(os.path.realpath(j),destfile)
      else:
        shutil.copy2(j,destfile)


  def create_surfaces(self):
    self.singleruns_names=dict()

    self.single_filenames={
            'neoin':'neo.in',
            'neo2in':'neo2.in',
            'neo2_exe':'create_surface.x',
            'in_file_axi': 'w7x-m24li.bc'} # Create Surf.py, #neo2.in
    self.single_filepaths=dict(
            neoin='/temp/gernot_k/Neo2/AG/Interface/Productive/2016_02_w7x-m24li/neo.in',
            neo2in= '/temp/gernot_k/Neo2/AG/Interface/Productive/2016_02_w7x-m24li/SURF1/neo2.in',
            neo2_exe='/temp/wakatobi/2020_Surface/new_scan/neo_2.x',
            in_file_axi= '/temp/gernot_k/Neo2/AG/Interface/Productive/2016_02_w7x-m24li/w7x-m24li.bc'
           )
    ##get_single_filenames

    s_vec, collpar_vec, Z_vec, T_vec = np.loadtxt(self._rundir+'surfaces.dat', unpack=True)

    k = 0
    for s in s_vec:
      #dir_name = surf_dir_name + str(k+1)
      dir_name = surf_dir_name + ('%.9e' %(s_vec[k])).replace('e', 'd').replace('d-', 'm')
      print(dir_name)
      self.singleruns_names[dir_name]={'boozer_s':s_vec[k],
                                      'conl_over_mfp':collpar_vec[k],
                                      'z_eff':Z_vec[k]}
      k = k + 1


  def createsinglefiles(self, overwrite=False, link=False):

    self.neo2nml=Neo2File(self.single_filepaths['neo2in'])

    for singlerun in self.singleruns_names:
      singlerunpath=os.path.join(self._rundir,singlerun)
      os.makedirs(singlerunpath,exist_ok=True)

      for i,j in self.single_filepaths.items():

        destfile=os.path.join(singlerunpath,os.path.basename(j))

        if os.path.isfile(destfile):
          if overwrite:
            print(destfile, ' is a file and is replaced')
            os.remove(destfile)
          else:
            print(destfile, ' is not updated') ### TODO Maybe check if it is the same file!!
            continue
        if i == 'neo2in':
          self.neo2nml.chval(**self.singleruns_names[singlerun])
          self.neo2nml.write(destfile)

        if link:
          os.symlink(os.path.realpath(j),destfile)
        else:
          shutil.copy2(j,destfile)
