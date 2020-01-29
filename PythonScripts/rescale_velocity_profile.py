#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Dec 12 2019

@author: Christopher Albert
"""

def rescale_velocity_profile(path:str, infilename:str, scaling_factor:float):
  """
  From a given neo2 hdf5-input file, create a new one with the velocity
  profile rescaled.

  input:
  ------
  path,str: string with path where the input file can be found.
  infilename: string with the name of the input file to use.
  scaling_factor: floating point number, scaling factor to apply to the
    velocity. Expected to be a single number applied to all radial
    positions.

  \todo Move to hdf5tools?
  """
  import h5py, hdf5tools
  import numpy as np

  ref = hdf5tools.get_hdf5file(path+infilename)

  vphi = ref['Vphi']
  vphi_rescaled = []

  for k in vphi:
    vphi_rescaled.append(k*scaling_factor)

  vphi_rescaled = np.array(vphi_rescaled)

  nameparts = infilename.rsplit('.', 1)

  vphi_namepart= '_vphi_{0}'.format(scaling_factor)

  if len(nameparts)==1:
    outfilename = nameparts[0]+vphi_namepart
  elif len(nameparts)==2:
    outfilename = nameparts[0]+vphi_namepart+'.'+nameparts[1]

  vphi_rs = hdf5tools.get_hdf5file_replace(outfilename)

  for k in ref.keys():
    ref.copy(source='/'+k, dest=vphi_rs, name='/' + k)

  dset = vphi_rs['Vphi']
  dset[...] = vphi_rescaled

  vphi_rs.close()
