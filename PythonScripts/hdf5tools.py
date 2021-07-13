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

def get_hdf5file_replace(filename: str):
  """Simple wrapper creating a new hdf5 file. Replaces existing file."""

  import h5py

  f = h5py.File(filename, 'w')

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
  values will be writen to a single numpy array.
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

def copy_hdf5_from_paths_to_single_file(paths: list, infilename: str, outfilename: str, ignore_hdf5_path: bool):
  """Combine files 'infilename' located at given 'paths' into a single file 'outfilename'.

  This function will collect hdf5 files with name 'infilename', at
  locations given via list 'paths', and write data into single file with
  name 'outfilename'.
  Data from a file is put into a group according to its path, unless
  'ignore_hdf5_path' is set to True, in which case data will be put into
  basename(normpath(path)). This is intended for combining results from
  different runs, to have only the flux surface label of the subfolders
  as a path inside the hdf5 file.

  input:
  ------
  paths: list of foldernames at which locations to look for input.
  infilename: name of the file to look for/from which to read data, i.e.
    it is assumed that all the files have the same name and are just in
    different paths.
  outfilename: name to use for file in which data is collected.
  ignore_hdf5_path: if False, data is copied from 'infilename' into
    'outfilename' at 'path', if true then location will be
    'basename(normpath(path))' instead.
  """
  from os.path import basename, join, normpath

  with get_hdf5file_replace(outfilename) as o:
    for path in paths:
      try:
        f = get_hdf5file(join(path, infilename))
      except OSError:
        print('No file ', infilename, ' found in ', path)
      else:
        if ignore_hdf5_path:
          f.copy(source='/', dest=o, name='/' + basename(normpath(path)))
        else:
          f.copy(source='/', dest=o, name='/' + path)
        f.close()

def copy_hdf5_from_subfolders_to_single_file(path, infilename: str, outfilename: str, ignore_hdf5_path: bool = False):
  """For 'infilename' in subfolders of 'path', join them into 'outfilename'.

  This will look for files named 'infilename' in subfolders of 'path'.
  There content will be joined into a single file named 'outfilename',
  placed into the current folder.
  The content of each in-file will be put into a group with the name of
  the subfolder.
  If a subfolder does not contain 'infilename', then the subfolder is
  ignored.

  input:
  ------
  path: either string with path to folders which contain data, or a list
    of such paths.
  infilename: name of the file(s) which to combine, i.e. which should be
    in subfolders of given path.
  outfilename: name of the file which should contain collected data.
  ignore_hdf5_path: passed to copy_from_paths_to_single_file, please
    look there for further information.
  """

  from itertools import chain
  from os import listdir
  from os.path import isfile, join

  folders = []

  # Make sure path is a list.
  if (type(path) == str):
    path = [path]

  for p in path:
    folders.append([join(p, f) for f in listdir(p) if not isfile(join(p, f))])

  folders = list(chain.from_iterable(folders))

  copy_hdf5_from_paths_to_single_file(folders, infilename, outfilename, ignore_hdf5_path)

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

  Note that there will be no ordering, it is implicitly assumed that
  elements are already ordered.

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

def collect_and_reshape(path: str, infilename: str, outfilename: str):
  """ Copy data to single file and reshape it.

  input:
  ------
  path: string with the path in which to find subfolders to check for
    input files.
  infilename: name of the input files.
  outfilename: name of the file into which to write the output, i.e. the
    collected and reshaped data.

  side effects:
  -------------
  Leaves a temporary file behind ('temp' + outfilename), from the copy
  process.
  """

  from os.path import join

  temp_filename = 'temp_' + outfilename
  copy_hdf5_from_subfolders_to_single_file(path, infilename, temp_filename)
  reshape_hdf5_file(join(path, temp_filename), join(path, outfilename))

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

def compare_hdf5_group_keys(reference_group, other_group, verbose: bool):
  """Check if two h5py groups contain the same subgroups/datasets.

  This function checks if two given h5py groups contain the same
  subgroups and datasets.
  The groups are considered to differ, if the lengths do not match, or
  if the reference does contain keys which are not in the other group.

  input
  ----------
  reference_group:
  other_group:

  return value
  ----------
  True, if the groups are the same, in the sense given above, false if not.

  side effects
  ----------
  There should be no side effects.

  limitations
  ----------
  - no indications how the files differ, i.e. which keys are additionaly
    in which file
  - No white/blacklisting. Required e.g. for final.h5 which adds group
    named after the folder in which the code runs.
  """
  import h5py

  lr = list(reference_group.keys())
  lo = list(other_group.keys())

  return_value = True

  if len(lr) != len(lo):
    # No longer considered to mean the file differ, new file might
    # include new data fields, not present in reference. The other way
    # round is already covered in the for-loop belop.
    # Thus just print a warning.
    # ~ return_value = False
    print('Warning: number of group keys differ in group ' + reference_group.name)

  for key in lr:
    if key in lo:
      if isinstance(reference_group[key], h5py.Group):
        return_value = return_value and compare_hdf5_group_keys(reference_group[key], other_group[key], verbose)
    else:
      return_value = False
      print("Key '" + key + "' only found in reference.")

  return return_value

def fill_list(list_):
  """Helper function: list or filename input to list output.

  Input of this helper function has the same requirements as for
  whitelist/blacklist parameter of compare_hdf5_group_data.
  The function then makes sure a list is returned, i.e. if the parameter
  contains an filename, it will read the content and return it as a
  list.
  If the argument passed is a string but the file could not be opened,
  an empty list is returned. A message about this is printed, but no
  exeption is thrown by this function.
  If the parameter contains already a list, it is simply returned.
  """
  if isinstance(list_, str):
    try:
      with open(list_) as f:
        l = f.readlines()
        return [x.strip() for x in l]
    except FileNotFoundError as e:
      print('File could not be opened, return empty list.')
      return list()
  elif isinstance(list_, list) and len(list_) > 0:
    return list_

def compare_hdf5_group_data(reference_group, other_group, delta_relative: float, whitelist, blacklist, verbose: bool):
  """Check if datsets of two h5py groups are equal.

  This function checks if two given h5py groups contain the same data,
  i.e. if the respective datasets are equal.
  Subgroups are checked recursively.
  Equality for a dataset of floats thereby means that the normalized
  values (to maximum value in dataset) differ less than a given delta.
  The difference is determined from the objects attribute 'accuracy' if
  it exists (considered to be an absolute accuracy) and from the
  parameter 'delta_relative' if the attribute does not exist or is
  either infinite or nan.
  For integers and strings it refers to direct equality.
  Either a whitelist of keys to compare or a blacklist of keys to not
  compare, can be given.

  input
  ----------
  reference_group, other_group: h5 (group) objects, the elements which
    should be compared.
  delta_relative: float, default value for the maximum error to use.
    This will only be used if the onject (1) is a float (2) has no
    attribute 'accuracy'. If the latter exists, then this is used as
    absolute accuracy.
  whitelist, blacklist: Both of these can be empty ([]), or one of them
    can either contain a string, or a list of strings.
    If it contains a string, this is assumed to be a filename which
    contains the entries for the whitelist/blacklist, one per line.
    If it contains a list, this is to be assumed to be the whitelist/
    blacklist.
    If both of them are not empty, this is considered an error.
  verbose: logical, if true, then the keys for which there are
    differences are printed, together with the relative error.

  return value
  ----------
  True, if the groups contain the same data, in the sense given above,
  false if not.

  side effects
  ----------
  There should be no side effects.

  limitations
  ----------
  Not sure how nan and inf are treated.
  At least encountering a nan should lead to a warning
  'RuntimeWarning: invalid value encountered in greater'.
  """
  import numpy
  import h5py
  from math import isnan, isinf

  lr = list(reference_group.keys())
  lo = list(other_group.keys())

  return_value = True

  if whitelist and blacklist:
    raise Exception()

  whitelist = fill_list(whitelist)
  blacklist = fill_list(blacklist)

  for key in lr:
    if (not whitelist or key in whitelist) and (not blacklist or key not in blacklist):
      if key in lo:
        if isinstance(reference_group[key], h5py.Dataset):
          if reference_group[key].dtype.kind == 'f':
            max_value_dataset = max(numpy.nditer(abs(numpy.array(reference_group[key]))))
            if max_value_dataset == 0:
              max_value_dataset = 1.0
            abs_differences = abs(numpy.subtract(numpy.array(reference_group[key]), numpy.array(other_group[key])))
            delta_abs = reference_group[key].attrs.get('accuracy', -1.0)
            if (delta_abs < 0) or isnan(delta_abs) or isinf(delta_abs):
              delta_abs = delta_relative * max_value_dataset

            if (abs_differences > delta_abs).any():
              return_value = False
              if (verbose):
                print('Difference in ' + key + ': ' + '{}'.format(float(max(numpy.nditer(abs_differences/max_value_dataset)))))
          elif reference_group[key].dtype.kind == 'i':
            if (numpy.array(reference_group[key]) != numpy.array(other_group[key])).any():
              return_value = False
              print('Difference in ' + key)
          elif reference_group[key].dtype.kind == 'S':
            if (numpy.array(reference_group[key]) != numpy.array(other_group[key])).any():
              return_value = False
              print('Difference in ' + key)

        else:
          return_value = return_value and compare_hdf5_group_data(reference_group[key], other_group[key], delta_relative, whitelist, blacklist, verbose)

  return return_value

def compare_hdf5_files(reference_filename: str, other_filename: str, delta_relative: float, whitelist, blacklist, verbose: bool):
  """Compare the content of two hdf5 files and return if they are equal or not.

  Compare the datasets of two hdf5 files to determine if they are equal
  within a certain accuracy.
  Two datasets are considered equal if |reference - other|/max(|reference|)
  is smaller than the given delta.
  The files are considered equal if this holds for all the datasets.

  input
  ----------
  reference_filename:
  other_filename:
  delta_relative:
  verbose: logical, if true, then the keys for which there are
    references are printed, together with the relative error.

  return value
  ----------
  List with to boolean values.
  The first one is true, if the files are the same, in the sense given
  above, false if not.
  The second one will be true of the keys of the two files match, i.e.
  the number and the names datasets are the same.

  side effects
  ----------
  There should be no side effects.

  limitations
  ----------
  - keys in reference need also to be in other (but not vice versa)
  """
  import numpy
  import h5py

  h5r = get_hdf5file(reference_filename)
  h5o = get_hdf5file(other_filename)

  lr = list(h5r.keys())
  lo = list(h5o.keys())

  keys_equal = compare_hdf5_group_keys(h5r, h5o, verbose)

  files_are_equal_to_delta = compare_hdf5_group_data(h5r, h5o, delta_relative, whitelist, blacklist, verbose)

  return [files_are_equal_to_delta, keys_equal]

def add_species_to_profile_file(infilename: str, outfilename: str, Zeff: float, Ztrace: float, mtrace: float):
  """Add to profile input file a species by scaling existing ion species.
  T_prof             needs to be changed ... done
  Vphi               ok
  boozer_s           ok
  dT_ov_ds_prof      needs to be changed ... done
  dn_ov_ds_prof      needs to be changed ... done
  isw_Vphi_loc       ok?
  kappa_prof         needs to be changed
  n_prof             needs to be changed ... done
  num_radial_pts     ok
  num_species        needs to be changed ... done
  rel_stages         needs to be changed ... done
  rho_pol            ok
  species_def        needs to be changed ... done
  species_tag        needs to be changed ... done
  species_tag_Vphi   ok?

  If Zeff = 1, then it is assumed that the density should be split
  evenly among the two ion species, as the term for Zeff /= 1 would
  fail. This might be usefull for having two hydrogen species, e.g. for
  testing.

  input
  ----------
  infilename:
  outfilename:
  Zeff: effective charge the plasma should have with the trace species.
  Ztrace: charge of the trace species in elementary charges.
  mtrace: mass of the trace species in [g]

  return value
  ----------
  No return value.

  side effects
  ----------
  New file is created.

  limitations
  ----------
  Assumes second species is ions.
  Assumes only a third species is added.
  """
  import math
  import numpy as np

  #~ Zeff = 1.3
  #~ Ztrace = 30
  #~ mtrace = 3.0527e-22 #[g]

  ELEMENTARY_CHARGE_SI = 1.60217662e-19
  SPEED_OF_LIGHT_SI = 2.99792458e8
  DENSITY_SI_TO_CGS = 1e-6
  ENERGY_SI_TO_CGS = 1e7
  EV_TO_SI = ELEMENTARY_CHARGE_SI
  EV_TO_CGS = EV_TO_SI * ENERGY_SI_TO_CGS
  CHARGE_SI_TO_CGS = 10*SPEED_OF_LIGHT_SI # Conversion factor is not speed of light, but 10c.

  if (Zeff == 1):
    factor_hydrogen = 0.5
    factor_trace = 0.5
  else:
    factor_hydrogen = (Ztrace - Zeff)/(Ztrace -1)
    factor_trace = (Zeff - 1)/(Ztrace*(Ztrace -1))

  no_change_needed = ['Vphi', 'boozer_s', 'isw_Vphi_loc', 'num_radial_pts', 'rho_pol', 'species_tag_Vphi']

  with get_hdf5file(infilename) as hin:
    with get_hdf5file_new(outfilename) as hout:
      for dname in no_change_needed:
        hout.create_dataset(dname, data=hin[dname])

      nsp = np.array(hin['num_species'])
      nsp[0] += 1
      hout.create_dataset('num_species', data=nsp)

      rs = np.array(hin['rel_stages'])
      rs[...] += 1
      hout.create_dataset('rel_stages', data=rs)

      t = np.array(hin['T_prof'])
      t.resize( (hout['num_species'][0], hout['num_radial_pts'][0]) )
      t[2, ...] = t[1, ...]
      hout.create_dataset('T_prof', data=t)

      dt = np.array(hin['dT_ov_ds_prof'])
      dt.resize( (hout['num_species'][0], hout['num_radial_pts'][0]) )
      dt[2, ...] = dt[1, ...]
      hout.create_dataset('dT_ov_ds_prof', data=dt)

      n = np.array(hin['n_prof'])
      n.resize( (hout['num_species'][0], hout['num_radial_pts'][0]) )
      n[1, ...] = factor_hydrogen * n[0, ...]
      n[2, ...] = factor_trace * n[0, ...]
      hout.create_dataset('n_prof', data=n)

      dn = np.array(hin['dn_ov_ds_prof'])
      dn.resize( (hout['num_species'][0], hout['num_radial_pts'][0]) )
      dn[1, ...] = factor_hydrogen * dn[0, ...]
      dn[2, ...] = factor_trace * dn[0, ...]
      hout.create_dataset('dn_ov_ds_prof', data=dn)

      st = np.array(hin['species_tag'])
      st.resize( (hout['num_species'][0], ) )
      st[2, ...] = hout['num_species'][0]
      hout.create_dataset('species_tag', data=st)

      sd = np.array(hin['species_def'])
      sd.resize((2, hout['num_species'][0], hout['num_radial_pts'][0]))
      # This reordering is necessary as not the added column is filled
      # with zeros.
      sd[1, 1, ...] = sd[1, 0, ...]
      sd[1, 0, ...] = sd[0, -1, ...]
      sd[0, -1, ...] = Ztrace
      sd[1, -1, ...] = mtrace
      hout.create_dataset('species_def', data=sd)

      Lambda = np.array([39.1 - 1.15*math.log10(x/DENSITY_SI_TO_CGS) + 2.3*math.log10(y/EV_TO_CGS*1e-3) for x, y in zip(n[0,...], t[0,...])])

      k = np.array(hin['kappa_prof'])
      k.resize( (hout['num_species'][0], hout['num_radial_pts'][0]) )
      # Value based on density and temperature, as the former changes,
      # needs to be recomputed for hydrogen, and computed for the trace.
      # As kappa seems to depent linearly on density and is independent
      # of charge, a simple rescaling might be enough?
      le = 3/(4*math.sqrt(math.pi)) * t[0,...]**2 / (n[0, ...] * (sd[0, 0, ...]*ELEMENTARY_CHARGE_SI*CHARGE_SI_TO_CGS)**4 * Lambda)
      li = 3/(4*math.sqrt(math.pi)) * t[1,...]**2 / (n[1, ...] * (sd[0, 1, ...]*ELEMENTARY_CHARGE_SI*CHARGE_SI_TO_CGS)**4 * Lambda)
      lt = 3/(4*math.sqrt(math.pi)) * t[2,...]**2 / (n[2, ...] * (sd[0, 2, ...]*ELEMENTARY_CHARGE_SI*CHARGE_SI_TO_CGS)**4 * Lambda)
      k[0, ...] = 2 / le
      k[1, ...] = 2 / li
      k[2, ...] = 2 / lt
      hout.create_dataset('kappa_prof', data=k)

def remove_species_from_profile_file(infilename: str, outfilename: str, index: int):
  """
  Function to remove a species from a profile input hdf5 file.

  This function is for removing a species from an input profile file.
  This might be usefull for testing, e.g. running just deuterium, or
  reducing the number of species to reduce memory requirements.

  input:
  ------
  infilename: file to read
  outfilename: name of the file to create
  index: zero based index of the species to remove

  output:
  -------
  No formal output. The result is written to a new hdf5 file with given
  name.
  """
  import numpy as np

  no_change_needed = ['Vphi', 'boozer_s', 'isw_Vphi_loc', 'num_radial_pts', 'rho_pol']
  special_treatment_needed = ['species_tag_Vphi']
  zero_dimension_species = ['T_prof', 'n_prof', 'kappa_prof', 'dn_ov_ds_prof', 'dT_ov_ds_prof']
  first_dimension_species = ['species_def']

  with get_hdf5file(infilename) as hin:
    with get_hdf5file_new(outfilename) as hout:
      for dname in no_change_needed:
        hout.create_dataset(dname, data=hin[dname])

      nspecies_tag_vphi = np.array(hin['species_tag_Vphi'])
      if index < nspecies_tag_vphi[0]:
        nspecies_tag_vphi[0] -= 1
      hout.create_dataset('species_tag_Vphi', data=nspecies_tag_vphi)

      nsp = np.array(hin['num_species'])
      nsp[0] -= 1
      hout.create_dataset('num_species', data=nsp)

      rs = np.array(hin['rel_stages'])
      rs[...] -= 1
      hout.create_dataset('rel_stages', data=rs)

      st = np.array(hin['species_tag'])
      st.resize( (hout['num_species'][0], ) )
      hout.create_dataset('species_tag', data=st)

      for dname in zero_dimension_species:
        t = np.array(hin[dname])[0:index, ...]
        t = np.append(t, np.array(hin[dname])[index+1:, ...], 0)
        hout.create_dataset(dname, data=t)

      for dname in first_dimension_species:
        t = np.array(hin[dname])[:, 0:index, ...]
        t = np.append(t, np.array(hin[dname])[:, index+1:, ...], 1)
        hout.create_dataset(dname, data=t)

def write_axisymmetric_quantities(infilename:str, outfilename:str):
  """ Write hardcoded set of output from hdf5-file to textfile.

  Was used for benchmarking with Jose Luis Velasco.
  """

  h = get_hdf5file(infilename)

  try:
    h['Gamma_NA_spec']
    header = "!               s             aiota          avnabpsi          Bref [G]                Er   MtOvR electrons   MtOvR deuterium D13  D23 D31 D32\n"
    er_available = True
  except KeyError:
    er_available = False
    header = "!               s             aiota          avnabpsi          Bref [G]\n"

  with open(outfilename, 'w') as outfile:
    outfile.write(header)
    for k in range(h['Bref'].size):
      if er_available:
        outfile.write("{: 16.10e} {: 16.10e} {: 16.10e} {: 16.10e} {: 16.10e} {: 16.10e} {: 16.10e} {: 16.10e} {: 16.10e} {: 16.10e} {: 16.10e} {: 16.10e} {: 16.10e} {: 16.10e}  {: 16.10e}\n".format(
          h['boozer_s'][k], h['aiota'][k], h['avnabpsi'][k], h['Bref'][k],
          h['Er'][k], h['MtOvR'][k][0], h['MtOvR'][k][1],
          h['D13_AX'][k][0], h['D13_AX'][k][3], h['D23_AX'][k][0],
          h['D23_AX'][k][3], h['D31_AX'][k][0], h['D31_AX'][k][3],
          h['D32_AX'][k][0], h['D32_AX'][k][3]))
      else:
        outfile.write("{: 16.10e} {: 16.10e} {: 16.10e} {: 16.10e} {: 16.10e} {: 16.10e} {: 16.10e} {: 16.10e} {: 16.10e} {: 16.10e} {: 16.10e} {: 16.10e}\n".format(
          h['boozer_s'][k], h['aiota'][k], h['avnabpsi'][k], h['Bref'][k],
          h['D13_AX'][k][0], h['D13_AX'][k][3], h['D23_AX'][k][0],
          h['D23_AX'][k][3], h['D31_AX'][k][0], h['D31_AX'][k][3],
          h['D32_AX'][k][0], h['D32_AX'][k][3]))

def write_nonaxisymmetric_quantities(infilename:str, outfilename:str):
  """ Write hardcoded set of output from hdf5-file to textfile.

  Was used for benchmarking with Jose Luis Velasco.
  """

  h = get_hdf5file(infilename)

  try:
    h['Gamma_NA_spec']
    header = "!               s    Gamma_NA elect     Gamma_NA deut  TphiNa electrons  TphiNa deuterium          D11ee_NA          D12ee_NA          D22ee_NA          D11ii_NA          D12ii_NA          D22ii_NA\n"
    gamma_available = True
  except KeyError:
    gamma_available = False
    header = "!               s          D11ee_NA          D12ee_NA          D22ee_NA          D11ii_NA          D12ii_NA          D22ii_NA\n"

  with open(outfilename, 'w') as outfile:
    outfile.write(header)
    for k in range(h['Bref'].size):
      if gamma_available:
        outfile.write("{: 16.10e} {: 16.10e} {: 16.10e} {: 16.10e} {: 16.10e} {: 16.10e} {: 16.10e} {: 16.10e} {: 16.10e} {: 16.10e} {: 16.10e}\n".format(
          h['boozer_s'][k], h['Gamma_NA_spec'][k][0], h['Gamma_NA_spec'][k][1], h['TphiNA_spec'][k][0], h['TphiNA_spec'][k][1],
          h['D11_NA'][k][0], h['D12_NA'][k][0], h['D22_NA'][k][0], h['D11_NA'][k][3], h['D12_NA'][k][3], h['D22_NA'][k][3]))
      else:
        outfile.write("{: 16.10e} {: 16.10e} {: 16.10e} {: 16.10e} {: 16.10e} {: 16.10e} {: 16.10e}\n".format(
          h['boozer_s'][k], h['D11_NA'][k][0], h['D12_NA'][k][0], h['D22_NA'][k][0], h['D11_NA'][k][3], h['D12_NA'][k][3], h['D22_NA'][k][3]))

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
