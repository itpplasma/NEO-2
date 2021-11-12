#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""\brief Wrapper script/function for running neo-2-par.

This is a wrapper function, that can also be used as script, to run
all four prop_reconstruct stages of a neo-2-par run

Example usage:
neo_2_par_wrapper.py -np nr_processors [last_stage_python] [name_executable]

Note that number of processors is here not an optional parameter.
Note that the number of processors is given with '-np' to remain
compatible with the condor submit file class, which at the moment
requires this for correct setting the number of processors.

For more information on the paramters check also the documentation of
the function 'neo_2_par_wrapper', but default values might differ.
"""


def neo_2_par_wrapper(nr_processors: int = 1, last_stage_python: bool = True, name_executable: str = "neo_2.x"):
  """\brief Wrapper function for running neo-2-par.

  input:
  ------
  nr_processors: optional integer, number of processors to use. Note that
    - mpi will be used if this is larger than three.
    - values of two will be ignored, and a single processor run is done
      instead.
    [1]
  last_stage_python: optional bool, if true then a python function from
    hdf5tools is used to collect the hdf5-files at the end. If false,
    also the last step is done with neo-2. [True]
  name_executable: str, name of the program to run. If num_processors>2
    then this is an argument to mpirun, for the parallel stages.
    ["neo_2.x"]
  """
  from os import system
  from pathlib import Path

  from scan_nonhdf5_tools import set_neo2in_reconstruction_for_folders
  from hdf5tools import prop_reconstruct_3

  multi_runcommand = 'mpirun -mca orte_tmpdir_base \"/tmp/\" -np {0} ./{1}'.format(nr_processors, name_executable)
  single_runcommand = "./" + name_executable

  folders = [Path('./')]

  runcommand = multi_runcommand
  if (nr_processors <= 2):
    runcommand = single_runcommand
    if (nr_processors == 2):
      print('WARNING: ignoring request for 2 processors, using only one!')

  set_neo2in_reconstruction_for_folders(folders, folder='./', backup=True, value=0)
  system(runcommand)
  set_neo2in_reconstruction_for_folders(folders, folder='./', backup=False, value=1)
  system(single_runcommand)
  set_neo2in_reconstruction_for_folders(folders, folder='./', backup=False, value=2)
  system(runcommand)
  if last_stage_python:
    prop_reconstruct_3(outfilename: str= 'final.h5'):
  else:
    set_neo2in_reconstruction_for_folders(folders, folder='./', backup=False, value=3)
    system(single_runcommand)


if __name__ == "__main__":
  import sys

  if (len(sys.argv) < 2):
    print("Usage:")
    print("./neo_2_par_wrapper.py -np nr_processors [last_stage_python] [name_executable]")
    print("  with radii in cm")
  else:
    nr_processors = 1
    if (len(sys.argv) >= 3):
      nr_processors = int(sys.argv[2])

    last_stage_python = True
    if (len(sys.argv) >= 4):
      last_stage_python = bool(sys.argv[3])

    name_executable = "neo_2.x"
    if (len(sys.argv) >= 5):
      name_executable = bool(sys.argv[4])

    neo_2_par_wrapper(nr_processors, last_stage_python, name_executable)
