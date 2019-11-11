#!/usr/bin/python

import os
import sys
import glob

basedir = os.path.dirname(os.path.realpath(__file__))
jobscript = os.path.join(basedir, 'neo2-single.job')
print(jobscript) 

subdirs = glob.glob('es_*')
subdirs.sort()

kjob = int(sys.argv[1])
njob = int(sys.argv[2])

endjob = kjob + njob 

for subdir in subdirs[kjob:endjob]:
  print('Submitting ' + subdir)
  os.chdir(subdir)
  os.system('sbatch ' + jobscript)
  os.chdir('..')

