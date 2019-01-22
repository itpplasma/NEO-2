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

def dim_zero(data, elements_to_keep: list, operate_on_last_dimension: bool):
  return data

def dim_one(data, elements_to_keep: list, operate_on_last_dimension: bool):
  import numpy as np
  return np.arry([data[elements_to_keep[x]] for x in range(len(elements_to_keep))])

def dim_two(data, elements_to_keep: list, operate_on_last_dimension: bool):
  import numpy as np
  if operate_on_last_dimension:
    return np.transpose(np.array([data[:, elements_to_keep[x]] for x in range(len(elements_to_keep))]))
  else:
    return np.array([data[elements_to_keep[x], :] for x in range(len(elements_to_keep))])

def dim_three(data, elements_to_keep: list, operate_on_last_dimension: bool):
  import numpy as np
  if operate_on_last_dimension:
    return np.moveaxis(np.array([data[:, :, elements_to_keep[x]] for x in range(len(elements_to_keep))]), 0, -1)
  else:
    return np.array([data[elements_to_keep[x], :, :] for x in range(len(elements_to_keep))])

def dim_four(data, elements_to_keep: list, operate_on_last_dimension: bool):
  import numpy as np
  if operate_on_last_dimension:
    return np.moveaxis(np.array([data[:, :, :, elements_to_keep[x]] for x in range(len(elements_to_keep))]), 0, -1)
  else:
    return np.array([data[elements_to_keep[x], :, :, :] for x in range(len(elements_to_keep))])

def dim_five(data, elements_to_keep: list, operate_on_last_dimension: bool):
  import numpy as np
  if operate_on_last_dimension:
    return np.moveaxis(np.array([data[:, :, :, :, elements_to_keep[x]] for x in range(len(elements_to_keep))]), 0, -1)
  else:
    return np.array([data[elements_to_keep[x], :, :, :, :] for x in range(len(elements_to_keep))])

def dim_six(data, elements_to_keep: list, operate_on_last_dimension: bool):
  import numpy as np
  if operate_on_last_dimension:
    return np.moveaxis(np.array([data[:, :, :, :, :,
      elements_to_keep[x]] for x in range(len(elements_to_keep))]), 0, -1)
  else:
    return np.array([data[elements_to_keep[x],
      :, :, :, :, :] for x in range(len(elements_to_keep))])

def dim_seven(data, elements_to_keep: list, operate_on_last_dimension: bool):
  import numpy as np
  if operate_on_last_dimension:
    return np.moveaxis(np.array([data[:, :, :, :, :,
      :, elements_to_keep[x]] for x in range(len(elements_to_keep))]), 0, -1)
  else:
    return np.array([data[elements_to_keep[x], :,
      :, :, :, :, :] for x in range(len(elements_to_keep))])

def dim_eigth(data, elements_to_keep: list, operate_on_last_dimension: bool):
  import numpy as np
  if operate_on_last_dimension:
    return np.moveaxis(np.array([data[:, :, :, :, :,
      :, :, elements_to_keep[x]] for x in range(len(elements_to_keep))]), 0, -1)
  else:
    return np.array([data[elements_to_keep[x], :, :, :,
      :, :, :, :, :] for x in range(len(elements_to_keep))])

def dim_nine(data, elements_to_keep: list, operate_on_last_dimension: bool):
  import numpy as np
  if operate_on_last_dimension:
    return np.moveaxis(np.array([data[:, :, :, :, :,
      :, :, :, elements_to_keep[x]] for x in range(len(elements_to_keep))]), 0, -1)
  else:
    return np.array([data[elements_to_keep[x], :, :, :,
      :, :, :, :, :] for x in range(len(elements_to_keep))])

def dim_ten(data, elements_to_keep: list, operate_on_last_dimension: bool):
  import numpy as np
  if operate_on_last_dimension:
    return np.moveaxis(np.array([data[:, :, :, :, :,
      :, :, :, :, elements_to_keep[x]] for x in range(len(elements_to_keep))]), 0, -1)
  else:
    return np.array([data[elements_to_keep[x], :, :, :, :,
      :, :, :, :, :] for x in range(len(elements_to_keep))])

def dim_eleven(data, elements_to_keep: list, operate_on_last_dimension: bool):
  import numpy as np
  if operate_on_last_dimension:
    return np.moveaxis(np.array([data[:, :, :, :, :,
      :, :, :, :, :,
      elements_to_keep[x]] for x in range(len(elements_to_keep))]), 0, -1)
  else:
    return np.array([data[elements_to_keep[x],
      :, :, :, :, :,
      :, :, :, :, :] for x in range(len(elements_to_keep))])

def dim_twelve(data, elements_to_keep: list, operate_on_last_dimension: bool):
  import numpy as np
  if operate_on_last_dimension:
    return np.moveaxis(np.array([data[:, :, :, :, :,
      :, :, :, :, :,
      :, elements_to_keep[x]] for x in range(len(elements_to_keep))]), 0, -1)
  else:
    return np.array([data[elements_to_keep[x], :,
      :, :, :, :, :,
      :, :, :, :, :] for x in range(len(elements_to_keep))])

def dim_thirteen(data, elements_to_keep: list, operate_on_last_dimension: bool):
  import numpy as np
  if operate_on_last_dimension:
    return np.moveaxis(np.array([data[:, :, :, :, :,
      :, :, :, :, :,
      :, :, elements_to_keep[x]] for x in range(len(elements_to_keep))]), 0, -1)
  else:
    return np.array([data[elements_to_keep[x], :, :,
      :, :, :, :, :,
      :, :, :, :, :] for x in range(len(elements_to_keep))])

def dim_fourteen(data, elements_to_keep: list, operate_on_last_dimension: bool):
  import numpy as np
  if operate_on_last_dimension:
    return np.moveaxis(np.array([data[:, :, :, :, :,
      :, :, :, :, :,
      :, :, :, elements_to_keep[x]] for x in range(len(elements_to_keep))]), 0, -1)
  else:
    return np.array([data[elements_to_keep[x], :, :, :,
      :, :, :, :, :,
      :, :, :, :, :] for x in range(len(elements_to_keep))])

def dim_fiveteen(data, elements_to_keep: list, operate_on_last_dimension: bool):
  import numpy as np
  if operate_on_last_dimension:
    return np.moveaxis(np.array([data[:, :, :, :, :,
      :, :, :, :, :,
      :, :, :, :, elements_to_keep[x]] for x in range(len(elements_to_keep))]), 0, -1)
  else:
    return np.array([data[elements_to_keep[x], :, :, :, :,
      :, :, :, :, :,
      :, :, :, :, :] for x in range(len(elements_to_keep))])

def remove_ends_from_hdf5_file(in_filename: str, out_filename: str,
      original_size: int, first_element_to_take: int,
      first_element_not_to_take: int, operate_on_last_dimension: bool):
  """Resize the arrays of an hdf5 file.

  This works like the version that takes a list instead of the index of
  the first element to take and the last element (not included).

  input
  ----------
  in_filename: File to read in and which content should be resized.
  out_filename: Name under which to store the resized file.
  original_size: The size of the (dimension of the) arrays which to
    resize. This is also used to determine if an array should be
    changed. For example if original_size=100 then 2x2 arrays (e.g.
    containing species data) will be left untouched.
  first_element_to_take: the zero based index of the first array item to
    keep.
  first_element_not_to_take: the zero based index of the array item up
    to which to keep. This index is not included, i.e. if given =10
    then only array items up to index 9 will be kept.
  operate_on_last_dimension: If this is true, then for multidimensional
    arrays, the last dimension will be resized. If false, then this
    function operates on the first dimension.

  return value
  ----------
  None, output is realized via side effect.

  side effects
  ----------
  Creates hdf5-file with given name 'out_filename', if first < last.
  """
  if first_element_to_take >= first_element_not_to_take:
    return

  resize_hdf5_file(in_filename, out_filename, original_size,
    range(first_element_to_take, first_element_not_to_take),
    operate_on_last_dimension)

def resize_hdf5_file(in_filename: str, out_filename: str, original_size: int, elements_to_keep: list, operate_on_last_dimension: bool):
  """Resize the arrays of an hdf5 file.

  input
  ----------
  in_filename: File to read in and which content should be resized.
  out_filename: Name under which to store the resized file.
  original_size: The size of the (dimension of the) arrays which to
    resize. This is also used to determine if an array should be
    changed. For example if original_size=100 then 2x2 arrays (e.g.
    containing species data) will be left untouched.
  elements_to_keep: list with the (zero based) indices of the elements,
    that should be keept.
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

  f = get_hdf5file(in_filename)
  o = get_hdf5file_new(out_filename)

  for key in list(f.keys()):
    size = list(f[key].shape)

    if operate_on_last_dimension:
      index_to_change = len(size)-1
    else:
      index_to_change = 0

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
    if (size[index_to_change] == original_size):
      dat = options[len(size)](f[key], elements_to_keep, operate_on_last_dimension)
    else:
      dat = f[key]

    o.create_dataset(key, data=dat)

  o['/num_radial_pts'][()] = np.array(len(elements_to_keep))

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
