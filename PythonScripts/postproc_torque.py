#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu, 2019-05-16

@author: Rico Buchholz

Based on shell/matlab script by Matyas Aradi.
"""

def data_process(folder: str, filename: str):
  from os.path import join
  from hdf5tools import get_hdf5file
  from numpy import array
  from math import pi
  import scipy.integrate as integrate

  # define constants
  e = 4.8032e-10 # elementary charge

  # load NEO-2 output
  with get_hdf5file(join(folder, filename)) as h5file:

    # allocate storage array
    aiota = []
    avbhat2 = []
    bcovar_phi = []
    boozer_psi_pr = []
    bcovar_tht = []
    boozer_s = []
    Bref = []
    MteOvR = []
    MtdOvR = []
    psi_pr_hat = []
    R0 = []
    TphiNA_ele = []
    TphiNA_io = []
    TphiNA_spec = []
    TphiNA_tot = []

    for k in list(h5file.keys()):
      if (not k.startswith('es_')):
        continue

      aiota.append(array(h5file[k]['aiota']).item())
      avbhat2.append(array(h5file[k]['avbhat2']).item())
      bcovar_phi.append(array(h5file[k]['bcovar_phi']).item())
      bcovar_tht.append(array(h5file[k]['bcovar_tht']).item())
      boozer_s.append(array(h5file[k]['boozer_s']).item())
      Bref.append(array(h5file[k]['Bref']).item())
      psi_pr_hat.append(array(h5file[k]['psi_pr_hat']).item())
      R0.append(array(h5file[k]['R0']).item())
      TphiNA_tot.append(array(h5file[k]['TphiNA_tot']).item())

      MtdOvR.append(h5file[k]['MtOvR'][1])
      MteOvR.append(h5file[k]['MtOvR'][0])
      TphiNA_spec.append(h5file[k]['TphiNA_spec'])

    Mte = array(MteOvR) * array(R0)
    Mtd = array(MtdOvR) * array(R0)
    avb2 = array(avbhat2) * array(Bref) * array(Bref)
    boozer_psi_pr = array(psi_pr_hat) * array(Bref);

    # Compute integral NTV torque
    TphiNA_ele = array(TphiNA_spec)[:,0]
    TphiNA_io = array(TphiNA_spec)[:,1]

    # compute integral torque
    TphiNA_int_tot = (4*pi*pi)*boozer_psi_pr[-1]*integrate.trapz(array(TphiNA_tot)*(array(aiota)*array(bcovar_tht)+array(bcovar_phi))/avb2, boozer_s)
    TphiNA_int_ele = (4*pi*pi)*boozer_psi_pr[-1]*integrate.trapz(TphiNA_ele*(array(aiota)*array(bcovar_tht)+array(bcovar_phi))/avb2, boozer_s)
    TphiNA_int_io = (4*pi*pi)*boozer_psi_pr[-1]*integrate.trapz(array(TphiNA_io)*(array(aiota)*array(bcovar_tht)+array(bcovar_phi))/avb2, boozer_s)

    return [boozer_s, TphiNA_tot, TphiNA_int_tot, TphiNA_int_ele, TphiNA_int_io, Mte, Mtd]

def plot_2spec_export(folder, vphifilename, outfilename, boozer_s, TphiNA_int_tot, TphiNA_int_ele, TphiNA_int_io):
  from os.path import join

  print('integral NTV torque = {:e} Nm'.format(1e-7*TphiNA_int_tot))
  print('integral NTV torque (ele) = {:e} Nm'.format(1e-7*TphiNA_int_ele))
  print('integral NTV torque (io) = {:e} Nm'.format(1e-7*TphiNA_int_io))

  vphiref = 0.0
  with open(join(folder, vphifilename)) as f:
    line = f.readline()
    parts = line.split('=')
    if (len(parts) > 1):
      vphiref = float(parts[1].split('!')[0])

  with open(join(folder, outfilename), 'w') as f:
    f.write('{:e} {:e} {:e} {:e}'.format(vphiref, 1e-7*TphiNA_int_tot, 1e-7*TphiNA_int_ele, 1e-7*TphiNA_int_io))

def export_2spec_Matyas(folder, h5filename, vphifilename, outfilename):
  [boozer_s, TphiNA_tot, TphiNA_int_tot, TphiNA_int_ele, TphiNA_int_io, Mte, Mtd] = data_process(folder, h5filename)
  plot_2spec_export(folder, vphifilename, outfilename, boozer_s, TphiNA_int_tot, TphiNA_int_ele, TphiNA_int_io)

def get_NTV_torque_int(folder: str, subfolder_pattern: str, filename: str):
  from os.path import join
  from pathlib import Path

  p = Path(folder)
  folders = list(p.glob(subfolder_pattern))

  content = []

  for d in folders:
    with open(join(folder, d.name, filename)) as f:
      content.append([])

      for numberstring in f.readline().split():
        content[-1].append(float(numberstring))

  content.sort()

  with open(join(folder, filename), 'w') as f:
    for numbers in content:
      f.write('{:+e} {:+e} {:+e} {:+e}\n'.format(numbers[0], numbers[1], numbers[2], numbers[3]))

def postproc_torque(folder: str, subfolder_pattern: str):
  from os.path import join
  from pathlib import Path
  from hdf5tools import copy_hdf5_from_subfolders_to_single_file

  infilename = 'neo2_multispecies_out.h5'

  p = Path(folder)
  folders = list(p.glob(subfolder_pattern))

  # for each folder wih a velocity shift ...
  for d in folders:
    print('Processing folder: ' + d.name)
    current_path_name = join(folder, d.name)

    # ... merge the neo2_multispecies output files ...
    copy_hdf5_from_subfolders_to_single_file(current_path_name, infilename, 'final_' + infilename)
    #print('- merging of hdf5 files done.')

    export_2spec_Matyas(current_path_name, 'final_neo2_multispecies_out.h5', 'vphiref.in', 'NTV_tot_test.dat')
    #print('- calculating NTV torque done.')

  get_NTV_torque_int(folder, subfolder_pattern, 'NTV_tot_test.dat')
