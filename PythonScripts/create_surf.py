#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Interactive script to create flux surface subfolders of neo-2-par run.
"""

if __name__ == "__main__":
  import numpy as np
  import shutil
  import subprocess

  surf_dir_name = "s"

  s_vec, collpar_vec, Z_vec, T_vec = np.loadtxt('surfaces.dat', unpack=True)

  print("Should I create the surface directories [y,n]?")
  ans = input().lower()
  print("Relativistic?")
  rel = input().lower()



  #ans = 'y'
  if ans == 'y':
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
      if rel == 'y':
        neo2 = neo2.replace('<T_e>',         ('%.9e' %(T_vec[k])).replace('e', 'd'))

      f = open(dir_name + "/neo2.in", "wt")
      f.write(neo2)
      f.close()
      k = k + 1

#print "Should I submit the jobs to condor [y,n]?"
#ans = raw_input().lower()
#if ans == 'y':
#  k = 0
#  for s in s_desired:
#    dir_name = surf_dir_name + str(k+1)
#    p = subprocess.Popen(["condor_submit_dag", "batch.dag"], cwd=dir_name)
#    p.wait()
#    k = k + 1


#print "Should I prepare jobfiles for DCluster [y,n]?"
#ans = raw_input().lower()
#if ans == 'y':
#  k = 0
#  for s in s_vec:
#    for j in range(0,3):
#      dir_name = surf_dir_name + str(k+1)
#      fin = open(dir_name + "/NEO2-"+str(j)+".job", "rt")
#      jobfile = fin.read()
#      jobfile = jobfile.replace('<jobname>', 'Neo2-' + dir_name)
#      fin.close()
#      fout = open(dir_name + "/NEO2-"+str(j)+".job", "wt")
#      fout.write(jobfile)
#      fout.close()
#    k = k + 1
