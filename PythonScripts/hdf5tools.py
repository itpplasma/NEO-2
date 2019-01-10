#!/usr/bin/env python3

"""Functions and definitions for handling hdf5 output of neo2.

This module contains functions and definitions for handling the hdf5
output of neo2. (Some functions/definitions might also be usefull in
general.)

Invokation as a script:

hdf5tools.py in_filename out_filename

e.g.

hdf5tools.py final_neo2_multispecies_out.h5 final_neo2_multispecies_out_over_boozer_s.h5

It will read in the file 'final_neo2_multispecies_out.h5', which is
expected to have at root only groups, which all have the same
dataelements. Each of the groups represents the result of a simulation
at a different value for a parameter.
The output file will have the same dataelements at root, as the
groups of the original file have, with the difference that they have now
an additional dimension, for the parameter.
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
    if (type(g[ksecond_level]) != type(g)):
      size = (len(keys), ) + g[ksecond_level].shape
      o.create_dataset('/' + ksecond_level, shape=size, dtype=g[ksecond_level].dtype)
      for attribute_name in list(g[ksecond_level].attrs.keys()):
        o['/' + ksecond_level].attrs[attribute_name] = g[ksecond_level].attrs[attribute_name]

  # Fill in the data.
  i = 0
  for ktosquash in list(f.keys()):
    g = f[ktosquash]
    for ksecond_level in list(o.keys()):
      dataset = o[ksecond_level]
      dataset[i] = g[ksecond_level]

    i = i + 1

def dim_zero(data, lower_bounds, upper_bounds):
  return data

def dim_one(data, lower_bounds, upper_bounds):
  return data[lower_bounds[0]:upper_bounds[0]]

def dim_two(data, lower_bounds, upper_bounds):
  return data[lower_bounds[0]:upper_bounds[0],
              lower_bounds[1]:upper_bounds[1]]

def dim_three(data, lower_bounds, upper_bounds):
  return data[lower_bounds[0]:upper_bounds[0],
              lower_bounds[1]:upper_bounds[1],
              lower_bounds[2]:upper_bounds[2]]

def dim_four(data, lower_bounds, upper_bounds):
  return data[lower_bounds[0]:upper_bounds[0],
              lower_bounds[1]:upper_bounds[1],
              lower_bounds[2]:upper_bounds[2],
              lower_bounds[3]:upper_bounds[3]]

def dim_five(data, lower_bounds, upper_bounds):
  return data[lower_bounds[0]:upper_bounds[0],
              lower_bounds[1]:upper_bounds[1],
              lower_bounds[2]:upper_bounds[2],
              lower_bounds[3]:upper_bounds[3],
              lower_bounds[4]:upper_bounds[4]]

def dim_six(data, lower_bounds, upper_bounds):
  return data[lower_bounds[0]:upper_bounds[0],
              lower_bounds[1]:upper_bounds[1],
              lower_bounds[2]:upper_bounds[2],
              lower_bounds[3]:upper_bounds[3],
              lower_bounds[4]:upper_bounds[4],
              lower_bounds[5]:upper_bounds[5]]

def dim_seven(data, lower_bounds, upper_bounds):
  return data[lower_bounds[0]:upper_bounds[0],
              lower_bounds[1]:upper_bounds[1],
              lower_bounds[2]:upper_bounds[2],
              lower_bounds[3]:upper_bounds[3],
              lower_bounds[4]:upper_bounds[4],
              lower_bounds[5]:upper_bounds[5],
              lower_bounds[6]:upper_bounds[6]]

def dim_eigth(data, lower_bounds, upper_bounds):
  return data[lower_bounds[0]:upper_bounds[0],
              lower_bounds[1]:upper_bounds[1],
              lower_bounds[2]:upper_bounds[2],
              lower_bounds[3]:upper_bounds[3],
              lower_bounds[4]:upper_bounds[4],
              lower_bounds[5]:upper_bounds[5],
              lower_bounds[6]:upper_bounds[6],
              lower_bounds[7]:upper_bounds[7]]

def dim_nine(data, lower_bounds, upper_bounds):
  return data[lower_bounds[0]:upper_bounds[0],
              lower_bounds[1]:upper_bounds[1],
              lower_bounds[2]:upper_bounds[2],
              lower_bounds[3]:upper_bounds[3],
              lower_bounds[4]:upper_bounds[4],
              lower_bounds[5]:upper_bounds[5],
              lower_bounds[6]:upper_bounds[6],
              lower_bounds[7]:upper_bounds[7],
              lower_bounds[8]:upper_bounds[8]]

def dim_ten(data, lower_bounds, upper_bounds):
  return data[lower_bounds[0]:upper_bounds[0],
              lower_bounds[1]:upper_bounds[1],
              lower_bounds[2]:upper_bounds[2],
              lower_bounds[3]:upper_bounds[3],
              lower_bounds[4]:upper_bounds[4],
              lower_bounds[5]:upper_bounds[5],
              lower_bounds[6]:upper_bounds[6],
              lower_bounds[7]:upper_bounds[7],
              lower_bounds[8]:upper_bounds[8],
              lower_bounds[9]:upper_bounds[9]]

def dim_eleven(data, lower_bounds, upper_bounds):
  return data[lower_bounds[0]:upper_bounds[0],
              lower_bounds[1]:upper_bounds[1],
              lower_bounds[2]:upper_bounds[2],
              lower_bounds[3]:upper_bounds[3],
              lower_bounds[4]:upper_bounds[4],
              lower_bounds[5]:upper_bounds[5],
              lower_bounds[6]:upper_bounds[6],
              lower_bounds[7]:upper_bounds[7],
              lower_bounds[8]:upper_bounds[8],
              lower_bounds[9]:upper_bounds[9],
              lower_bounds[10]:upper_bounds[10]]

def dim_twelve(data, lower_bounds, upper_bounds):
  return data[lower_bounds[0]:upper_bounds[0],
              lower_bounds[1]:upper_bounds[1],
              lower_bounds[2]:upper_bounds[2],
              lower_bounds[3]:upper_bounds[3],
              lower_bounds[4]:upper_bounds[4],
              lower_bounds[5]:upper_bounds[5],
              lower_bounds[6]:upper_bounds[6],
              lower_bounds[7]:upper_bounds[7],
              lower_bounds[8]:upper_bounds[8],
              lower_bounds[9]:upper_bounds[9],
              lower_bounds[10]:upper_bounds[10],
              lower_bounds[11]:upper_bounds[11]]

def dim_thirteen(data, lower_bounds, upper_bounds):
  return data[lower_bounds[0]:upper_bounds[0],
              lower_bounds[1]:upper_bounds[1],
              lower_bounds[2]:upper_bounds[2],
              lower_bounds[3]:upper_bounds[3],
              lower_bounds[4]:upper_bounds[4],
              lower_bounds[5]:upper_bounds[5],
              lower_bounds[6]:upper_bounds[6],
              lower_bounds[7]:upper_bounds[7],
              lower_bounds[8]:upper_bounds[8],
              lower_bounds[9]:upper_bounds[9],
              lower_bounds[10]:upper_bounds[10],
              lower_bounds[11]:upper_bounds[11],
              lower_bounds[12]:upper_bounds[12]]

def dim_fourteen(data, lower_bounds, upper_bounds):
  return data[lower_bounds[0]:upper_bounds[0],
              lower_bounds[1]:upper_bounds[1],
              lower_bounds[2]:upper_bounds[2],
              lower_bounds[3]:upper_bounds[3],
              lower_bounds[4]:upper_bounds[4],
              lower_bounds[5]:upper_bounds[5],
              lower_bounds[6]:upper_bounds[6],
              lower_bounds[7]:upper_bounds[7],
              lower_bounds[8]:upper_bounds[8],
              lower_bounds[9]:upper_bounds[9],
              lower_bounds[10]:upper_bounds[10],
              lower_bounds[11]:upper_bounds[11],
              lower_bounds[12]:upper_bounds[12],
              lower_bounds[13]:upper_bounds[13]]

def dim_fiveteen(data, lower_bounds, upper_bounds):
  return data[lower_bounds[0]:upper_bounds[0],
              lower_bounds[1]:upper_bounds[1],
              lower_bounds[2]:upper_bounds[2],
              lower_bounds[3]:upper_bounds[3],
              lower_bounds[4]:upper_bounds[4],
              lower_bounds[5]:upper_bounds[5],
              lower_bounds[6]:upper_bounds[6],
              lower_bounds[7]:upper_bounds[7],
              lower_bounds[8]:upper_bounds[8],
              lower_bounds[9]:upper_bounds[9],
              lower_bounds[10]:upper_bounds[10],
              lower_bounds[11]:upper_bounds[11],
              lower_bounds[12]:upper_bounds[12],
              lower_bounds[13]:upper_bounds[13],
              lower_bounds[14]:upper_bounds[14]]

def resize_hdf5_file(in_filename: str, out_filename: str, original_size: int, first: int, last: int, operate_on_last_dimension: bool):
  """Resize the arrays of an hdf5 file.

  input
  ----------
  in_filename: File to read in and which content should be resized.
  out_filename: Name under which to store the resized file.
  original_size: The size of the (dimension of the) arrays which to
    resize. This is also used to determine if an array should be
    changed. For example if original_size=100 then 2x2 arrays (e.g.
    containing species data) will be left untouched.
  first: the zero based index of the first array item to keep.
  last: the zero based index of the array item up to which to keep. This
    index is not included, i.e. if given last=10 then only array items
    up to index 9 will be kept.
  operate_on_last_dimension: If this is true, then for multidimensional
    arrays, the last dimension will be resized. If false, then this
    function operates on the first dimension.

  return value
  ----------
  None, output is realized via side effect.

  side effects
  ----------
  Creates hdf5-file with given name 'out_filename', if first < last.

  limitations
  ----------
  So far it is assumed that the key with the number of points for the
  resized dimension is '/num_radial_pts'.
  """

  import h5py
  import numpy as np

  # If the 'list' of indices is empty, do nothing.
  if first >= last:
    return

  f = get_hdf5file(in_filename)
  o = get_hdf5file_new(out_filename)

  for key in list(f.keys()):
    size = list(f[key].shape)

    lower_bounds = [0 for x in range(len(size))]
    upper_bounds = [x for x in size]

    if operate_on_last_dimension:
      index_to_change = len(size)-1
    else:
      index_to_change = 0

    if (original_size == size[index_to_change]):
      size[index_to_change] = last - first
      lower_bounds[index_to_change] = first
      upper_bounds[index_to_change] = last

    #~ dat = o.create_dataset(key, shape=tuple(size), dtype=f[key].dtype)

    # Fill the dataset with values. This codes uses a suggestion from
    #https://stackoverflow.com/questions/11479816/what-is-the-python-equivalent-for-a-case-switch-statement/11479840#11479840
    # to recreate a switch/select statement in python.

    # map the inputs to the function blocks
    options = {0 : dim_zero,
               1 : dim_one,
               2 : dim_two,
               3 : dim_three,
               4 : dim_four,
               5 : dim_five,
               6 : dim_six,
               7 : dim_seven,
               8 : dim_eigth,
               9 : dim_nine,
               10: dim_ten,
               11: dim_eleven,
               12: dim_twelve,
               13: dim_thirteen,
               14: dim_fourteen,
               15: dim_fiveteen,
    }

    # Call the corresponding function, which corresponds to the body of
    # the 'case'.
    #~ dat = options[len(size)](f[key], lower_bounds, upper_bounds)
    dat = options[len(size)](f[key], lower_bounds, upper_bounds)

    o.create_dataset(key, data=dat)

  o['/num_radial_pts'][()] = np.array([last - first])

if __name__ == "__main__":

  import h5py
  import sys

  # Check input arguments for input filename.
  if (len(sys.argv) >= 2):
    in_filename = sys.argv[1]
  else:
    in_filename = 'final_neo2_multispecies_out.h5'

  if (len(sys.argv) >= 3):
    out_filename = sys.argv[2]
  else:
    out_filename = 'final_neo2_multispecies_out_over_boozer_s.h5'

  reshape_hdf5_file(in_filename, out_filename)
