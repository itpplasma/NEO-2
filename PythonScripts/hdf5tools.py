#!/usr/bin/env python3

"""Functions and definitions for handling hdf5 output of neo2.

This module contains functions and definitions for handling the hdf5
output of neo2. (Some functions/definitions might also be usefull in
general.)

Invokation as a script:

hdf5tools.py fileforyaxis dataelementy fileforxaxis dataelementx

e.g.

hdf5tools.py fulltransp.h5 gamma_out neo2_config.h5 collision/conl_over_mfp

It will read from subfolders the element 'dataelementy' of the hdf5 file
'fileforyaxis' and analog for 'dataelementx' and 'fileforxaxis'.
It is assumed that 'dataelementy' is a matrix with at least three
dimensions in each direction as plots for these three times three
elements are done.
"""

# From h5py documentation:
#   Groups work like dictionaries, and datasets work like NumPy arrays

def get_hdf5file(filename: str):
  """Simple wrapper for reading a hdf5 file."""

  import h5py

  f = h5py.File(filename, 'r')

  return f

def get_hdf5file_to_append(filename: str):
  """Simple wrapper for appending to a hdf5 file."""

  import h5py

  f = h5py.File(filename, 'a')

  return f

def get_hdf5file_new(filename: str):
  """Simple wrapper creating a new hdf5 file. Will fail if file exists."""

  import h5py

  f = h5py.File(filename, 'w-')

  return f

def get_hdf5data(hdf5file, dataname: str):
  """Wrapper for getting data from hdf5 file object returned by get_hdf5file."""

  import h5py

  # Get the data
  data = hdf5file[dataname]

  return data

def get_hdf5data_from_subfolders(path: str, filename: str, dataname: str):
  """Get dataset from each 'filename' in subfolders of 'path'.

  This will check the given 'path' for folders. In each folder it will
  look for the hdf5-file 'filename' and read 'dataname' from it. The
  values will be writeen to a single numpy array.
  """

  import numpy as np
  from os import listdir
  from os.path import isfile, join

  values = []

  folders = [f for f in listdir(path) if not isfile(join(path, f))]

  for foldername in folders:
    #~ print(foldername)
    f = get_hdf5file(join(foldername, filename))
    d = np.array(get_hdf5data(f, dataname))
    values.append(d)

  return np.array(values)

def copy_hdf5_from_subfolders_to_single_file(path: str, infilename: str, outfilename: str):
  """For 'infilename' in subfolders of 'path', join them into 'outfilename'.

  This will look for files named 'infilename' in subfolders of 'path'.
  There content will be joined into a single file named 'outfilename'.
  The content of each in-file will be put into a group with the name of
  the subfolder.
  """

  from os import listdir
  from os.path import isfile, join

  values = []

  folders = [f for f in listdir(path) if not isfile(join(path, f))]

  try:
    o = get_hdf5file_new(outfilename)
  except OSError as error:
    print('OSError: {0}'.format(error))
    raise error
  else:
    for foldername in folders:
      try:
        f = get_hdf5file(join(foldername, infilename))
      except OSError:
        print('No file ', infilename,' found in ', foldername)
      else:
        f.copy(source='/', dest=o, name='/' + foldername)
        f.close()

def sort_x_y_pairs_according_to_x(x: list, y: list):
  """Sort pair of lists according to the first.

  This function takes to two lists, and will sort both lists according
  to the entries of the first list. You can think of this as either
  the same permutation (which sorts the entries in the first list) is
  applied to both, or you can say you create a list of pairs, sort
  according to the x-part of the pairs, and then make two seperate lists
  again.
  """

  intermediate = list(zip(x,y))

  intermediate = sorted(intermediate)

  x2, y2 = zip(*intermediate)

  return [np.array(x2), np.array(y2)]

def reshape_hdf5_file(in_filename: str, out_filename: str):
  """Reshape hdf5 file from x-times the 'same' group to group of arrays.

  This function assumes that 'filename' is the name of a hdf5 file,
  which contains mulitiple groups at the root, which all have the same
  structure, i.e. the names and forms of subgroups are the same.
  The output will be a file with one group with the same structure as
  those of a subgroup from the original file, but where each data
  element has now one dimension more, which contains the original
  top-level group.
  Example (the existence of subsubgoups is not implemented):
  /a1
    /b 1.0
    /c 1.0
       0.9
  /a2
    /b 1.1
    /c 2.0
       0.8
  /a3
    /b 1.2
    /c 3.0
       0.7
  will become
  /b
    1.0
    1.1
    1.2
  /c
    1.0 0.9
    2.0 0.8
    3.0 0.7
  """

  import h5py
  f = get_hdf5file(in_filename)

  o = get_hdf5file_to_append(out_filename)

  # Create basic file layout, including the sizes.
  keys = list(f.keys())
  g = f[keys[1]]
  for ksecond_level in list(g.keys()):
    size = (len(keys), ) + g[ksecond_level].shape
    o.create_dataset('/' + ksecond_level, shape=size, dtype=g[ksecond_level].dtype)

  # Fill in the data.
  i = 0
  for ktosquash in list(f.keys()):
    g = f[ktosquash]
    for ksecond_level in list(g.keys()):
      dataset = o[ksecond_level]
      dataset[i] = g[ksecond_level]

    i = i + 1

# Example, for how to call this as a script:
#   hdf5tools.py fulltransp.h5 gamma_out neo2_config.h5 collision/conl_over_mfp
#
# This would be the same, as calling this as a script without arguments.
# The first parameter is the file from which to take the second
# parameter (i.e. this is the name of the variable which to plot).
# The third and fourth parameter fulfill similar roles, just for the
# values of the x-axes.
# This means as default the field gamma_out is taken from the file
# fulltransp.h5 and ploted versus conl_over_mfp from the group collision
# in the file neo2_config.h5.
# What can not be seen from the example, is that the field gamma_out
# is a 3x3 array, and each component will be plotted and saved to a
# different file.
if __name__ == "__main__":

  import h5py
  import matplotlib.pyplot as plt
  import numpy as np
  import sys

  # Check input arguments for input filename.
  if (len(sys.argv) >= 2):
    filename = sys.argv[1]
  else:
    filename = 'fulltransp.h5'

  if (len(sys.argv) >= 3):
    fieldname = sys.argv[2]
  else:
    fieldname = 'gamma_out'

  if (len(sys.argv) >= 4):
    filename_x = sys.argv[3]
  else:
    filename_x = 'neo2_config.h5'

  if (len(sys.argv) >= 5):
    fieldname_x = sys.argv[4]
  else:
    fieldname_x = 'collision/conl_over_mfp'

  d = get_hdf5data_from_subfolders('./', filename, fieldname)
  x = get_hdf5data_from_subfolders('./', filename_x, fieldname_x)

  [x, d] = sort_x_y_pairs_according_to_x(x,d)

  plt.figure()
  plt.plot(x, d[:,1-1,1-1])
  plt.xscale('log')
  plt.yscale('log')
  plt.xlabel(fieldname_x)
  plt.ylabel(fieldname + ' 11')
  plt.savefig('gamma_11')

  plt.figure()
  plt.plot(x, d[:,1-1,2-1])
  plt.xscale('log')
  plt.yscale('log')
  plt.xlabel(fieldname_x)
  plt.ylabel(fieldname + ' 12')
  plt.savefig('gamma_12')

  plt.figure()
  plt.plot(x, d[:,1-1,3-1])
  plt.xscale('log')
  plt.yscale('symlog', linthreshx=1e-9)
  plt.xlabel(fieldname_x)
  plt.ylabel(fieldname + ' 13')
  plt.savefig('gamma_13')

  plt.figure()
  plt.plot(x, d[:,2-1,1-1])
  plt.xscale('log')
  plt.yscale('log')
  plt.xlabel(fieldname_x)
  plt.ylabel(fieldname + ' 21')
  plt.savefig('gamma_21')

  plt.figure()
  plt.plot(x, d[:,2-1,2-1])
  plt.xscale('log')
  plt.yscale('log')
  plt.xlabel(fieldname_x)
  plt.ylabel(fieldname + ' 22')
  plt.savefig('gamma_22')

  plt.figure()
  plt.plot(x, d[:,2-1,3-1])
  plt.xscale('log')
  plt.yscale('symlog', linthreshx=1e-9)
  plt.xlabel(fieldname_x)
  plt.ylabel(fieldname + ' 23')
  plt.savefig('gamma_23')

  plt.figure()
  plt.plot(x, d[:,3-1,1-1])
  plt.xscale('log')
  plt.yscale('symlog', linthreshx=1e-9)
  plt.xlabel(fieldname_x)
  plt.ylabel(fieldname + ' 31')
  plt.savefig('gamma_31')

  plt.figure()
  plt.plot(x, d[:,3-1,2-1])
  plt.xscale('log')
  plt.yscale('symlog', linthreshx=1e-9)
  plt.xlabel(fieldname_x)
  plt.ylabel(fieldname + ' 32')
  plt.savefig('gamma_32')

  plt.figure()
  plt.plot(x, d[:,3-1,3-1])
  plt.xscale('log')
  plt.yscale('symlog', linthreshx=1e-9)
  plt.xlabel(fieldname_x)
  plt.ylabel(fieldname + ' 33')
  plt.savefig('gamma_33')
