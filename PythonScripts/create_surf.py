#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to create flux surface subfolders of neo-2-par run.
"""

def create_surfaces(surface_dat_file:str = 'surfaces.dat'):
  import numpy as np
  import shutil

  surf_dir_name = "s"

  s_vec, collpar_vec, Z_vec, T_vec = np.loadtxt(surface_dat_file, unpack=True)

  if isinstance(s_vec, float): # if only one surface is given
    s_vec = np.array([s_vec])
    collpar_vec = np.array([collpar_vec])
    Z_vec = np.array([Z_vec])
    T_vec = np.array([T_vec])

  k = 0
  for s in s_vec:
    #dir_name = surf_dir_name + str(k+1)
    dir_name = surf_dir_name + ('%.9e' %(s_vec[k])).replace('e', 'd').replace('d-', 'm')
    print(dir_name)
    shutil.rmtree(dir_name, ignore_errors=True)
    shutil.copytree('TEMPLATE_DIR', dir_name, symlinks=True)
    f = open(dir_name + "/neo2.in", "rt")
    neo2 = f.read()
    f.close()

    neo2 = neo2.replace('<conl_over_mfp>', ('%.9e' %(collpar_vec[k])).replace('e', 'd'))
    neo2 = neo2.replace('<boozer_s>',      ('%.9e' %(s_vec[k])).replace('e', 'd'))
    neo2 = neo2.replace('<z_eff>',         ('%.9e' %(Z_vec[k])).replace('e', 'd'))
    neo2 = neo2.replace('<T_e>',         ('%.9e' %(T_vec[k])).replace('e', 'd'))

    f = open(dir_name + "/neo2.in", "wt")
    f.write(neo2)
    f.close()
    k = k + 1


if __name__ == "__main__":
  create_surfaces()
