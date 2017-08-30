#! /usr/bin/env python
import os
import f90nml
import shutil

for dirpath, dirnames, files in os.walk('./'):
    precom=False
    if dirpath=="./":
        if "f90nml" not in dirnames:
            print("Fortran Python Modul f90nml must be add to init_rund Directory")
            break
        if "hosts" not in files:
            print("hosts File in init_run Directory needed")
            break
        continue
        
    if "neo2_multispecies_out.h5" in files: # Look for finished Files
        print(dirpath + " is finished")
        continue
    if "neo2.in" in files:
        neo2=f90nml.read("neo2.in")
        
 
        
        if "precom_collop.h5" in files:
            neo2['collision']['lsw_write_precom']=False
            neo2['collision']['lsw_read_precom']=True
            neo2.write("neo2.in",force=True)
            print("neo2.in changed from Write to Read in " + dirpath)
            precom=True
        else:
            print("precom_collop.h5 does not exist in " + dirpath )
        
        if "hosts" not in files:
            
            shutil.copy('hosts',dirpath+'/hosts')
            print("Added hostfile to "+ dirpath)
        else:
            if precom:
                print(dirpath + " is ready to start with existing precom_collop")