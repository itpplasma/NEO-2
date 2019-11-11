#!/usr/bin/python

import os
import glob

proc_id = int(os.environ['SLURM_PROCID'])
print(proc_id)

subdirs = glob.glob('es_*')
subdirs.sort()
subdir = subdirs[proc_id]
print(subdir)

os.chdir(subdir)
os.system('module load mkl impi && LD_LIBRARY_PATH=$MKLROOT/lib/intel64:$LD_LIBRARY_PATH mpiexec -np 2 ./neo_2.x')

