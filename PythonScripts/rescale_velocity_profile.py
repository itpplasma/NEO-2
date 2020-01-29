#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Dec 12 2019

@author: Christopher Albert
"""

import h5py, hdf5tools
import numpy as np

ref = hdf5tools.get_hdf5file('multi_spec_aug32169_t4.1500.in')

vphi = ref['Vphi']
vphi_lower = []
vphi_larger = []

for k in vphi:
  vphi_lower.append(k*0.9)
  vphi_larger.append(k*1.1)

vphi_lower = np.array(vphi_lower)
vphi_larger = np.array(vphi_larger)

vphi_lo = hdf5tools.get_hdf5file_replace('multi_spec_aug32169_t4.1500_vphi_09.in')
vphi_la = hdf5tools.get_hdf5file_replace('multi_spec_aug32169_t4.1500_vphi_11.in')

for k in ref.keys():
  ref.copy(source='/'+k, dest=vphi_lo, name='/' + k)
  ref.copy(source='/'+k, dest=vphi_la, name='/' + k)

dset = vphi_lo['Vphi']
dset[...] = vphi_lower
dset = vphi_la['Vphi']
dset[...] = vphi_larger

vphi_lo.close()
vphi_la.close()
