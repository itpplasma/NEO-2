#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu, 2019-05-16

@author: Rico Buchholz

Based on shell/matlab script by Matyas Aradi.
"""

def data_process(folder: str, filename: str):
  """Process output of a scan over radius to get the torque.

  This function will process the output data (specifically, some fields
  from 'final_neo2_multispecies_out.h5', i.e. the combined data from the
  single runs, to compute the integral torque in cgs units.

  Note: there are some additional output quantities. Not sure what to do
  with these.

  input:
  ------
  folder: string, location of the radial scan.
  filename: string, name of the hdf5 file (not path), that contains the
    merged data of the radial scan.
  """
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

def plot_2spec_export(folder, vphifilename, boozer_s, TphiNA_int_tot, TphiNA_int_ele, TphiNA_int_io):
  """Combine the torque data with vphiref and convert to SI units.

  This function will read the vphiref value from file. The calculated
  torques are printed in SI units. Torques (in SI units) and vphiref are
  returned.

  input:
  ------
  folder: string, path where the vphifile is.
  vphifilename: string, name of the file, that contains the value of
    vphiref.
  boozer_s: list with s values of radial scan.
  TphiNA_int_tot: floating point number, with the value of the
    integrated total torque.
  TphiNA_int_ele: floating point number, with the value of the
    integrated electron torque.
  TphiNA_int_io: floating point number, with the value of the
    integrated ion torque.
  """
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

  return [vphiref, 1e-7*TphiNA_int_tot, 1e-7*TphiNA_int_ele, 1e-7*TphiNA_int_io]

def export_2spec_Matyas(folder, h5filename, vphifilename):
  """Process data of radial scan and return vphiref and torque (in SI units).

  Process data of radial scan to calculate torque in SI units and return
  them together with vphiref.

  input:
  ------
  folder: string, path where the data of the radial scan is.
  h5filename: string, name of the h5file (not path), that contains the
    (merged) data of the radial scan.
  vphifilename: string, name of the file (not path), that contains the
    vphiref value of the radial scan.
  """
  [boozer_s, TphiNA_tot, TphiNA_int_tot, TphiNA_int_ele, TphiNA_int_io, Mte, Mtd] = data_process(folder, h5filename)
  return plot_2spec_export(folder, vphifilename, boozer_s, TphiNA_int_tot, TphiNA_int_ele, TphiNA_int_io)

def write_torque_data(folder: str, outfilename: str, torque_data):
  """Write data ordered to a file.

  This function assumes that torque_data is a list of lists with four
  floating point entries.
  The data will be sorted according to the first entries of the inner
  list and then written to a file 'outfilename' located in 'folder'.
  Each inner list will be put into one line.

  input:
  ------
  folder: string, path where to put the output file.
  outfilename: string, name of the output file.
  torque_data: list of vphi and torque values, realized as list of
    lists.
  """
  from os.path import join
  from pathlib import Path

  # This should sort the list according to the first entries of the inner list.
  torque_data.sort()

  with open(join(folder, outfilename), 'w') as f:
    for numbers in torque_data:
      f.write('{:+e} {:+e} {:+e} {:+e}\n'.format(numbers[0], numbers[1], numbers[2], numbers[3]))

def postproc_torque(folder: str, subfolder_pattern: str, hdf5infilename: str, vphiinfilename: str, torqueoutfilename: str):
  """Process data of scans over vphiref, calculate the torque and write it to file.

  This function will process data from a scan over vphiref. First the
  hdf5 output from the radial scan is merged into a single file per
  vphiref value. From these the torque is calculated in SI units. The
  torque and vphiref is written to file.

  Example, if the subfolders are in the current folder:
  postproc_torque('./', 'n2_vshift*', 'neo2_multispecies_out.h5', 'vphiref.in', 'NTV_torque_int_all.dat')

  input:
  ------
  folder: string, working directory where the scan over vphiref is
    located. Use './' if the current directory is the working directory.
  subfolder_pattern: string, pattern to use for the subfolders, which
    contain the individual runs of the vphiref scan. Example:
    'n2_vshift+*'.
  hdf5infilename: string, name of the hdf5 file, which contains the data
    of a single vphiref and s run.
    Attention: this is only the name not the path.
  vphiinfilename: string, name of the file, that contains the value of
    vphiref in rad/s.
    Attention: this is only the name not the path.
  torqueoutfilename: string, name of the file, where to store the
    torque.
    Attention: this is only the name not the path.
  """
  from os.path import join
  from pathlib import Path
  from hdf5tools import copy_hdf5_from_subfolders_to_single_file

  p = Path(folder)
  folders = list(p.glob(subfolder_pattern))
  # Note: torque_data might not be ordered according to vphi.
  torque_data = []

  # for each folder wih a velocity shift ...
  for d in folders:
    print('Processing folder: ' + d.name)
    current_path_name = join(folder, d.name)

    # ... merge the neo2_multispecies output files ...
    copy_hdf5_from_subfolders_to_single_file(current_path_name, hdf5infilename, 'final_' + hdf5infilename)
    #print('- merging of hdf5 files done.')

    torque_data.append(export_2spec_Matyas(current_path_name, 'final_' + hdf5infilename, vphiinfilename))
    #print('- calculating NTV torque done.')

  write_torque_data(folder, torqueoutfilename, torque_data)
