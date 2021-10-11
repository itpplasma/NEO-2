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

def copy_hdf5_from_paths_to_single_file(paths: list, infilename: str, outfilename: str, only_base_path: bool, source_base_path: bool):
  """Combine files 'infilename' located at given 'paths' into a single file 'outfilename'.

  This function will collect hdf5 files with name 'infilename', at
  locations given via list 'paths', and write data into single file with
  name 'outfilename'.
  Data from a file is put into a group according to its path, unless
  'only_base_path' is set to True, in which case data will be put into
  basename(normpath(path)). This is intended for combining results from
  different runs, to have only the flux surface label of the subfolders
  as a path inside the hdf5 file.
  Note that data can not simply be copied from root to root, as this
  would lead to an exception, as the group already exists.
  The user can decide to copy everything, i.e. from root, or just those
  data that is in group '/'+'basename(normpath(path))'. The latter is
  thought for collecting files like 'final.h5' as created by neo-2-par
  at the end of a reconstruction run.

  input:
  ------
  paths: list of foldernames at which locations to look for input.
  infilename: name of the file to look for/from which to read data, i.e.
    it is assumed that all the files have the same name and are just in
    different paths.
  outfilename: name to use for file in which data is collected.
  only_base_path: if False, data is copied from root ('/') of
    'infilename' into a group named 'path' in 'outfilename', if True
    then groupname will be 'basename(normpath(path))' instead.
  source_base_path: if False, copy everything, if True then only copy
    from group '/'+'basename(normpath(path))'.
  """
  from os.path import basename, join, normpath

  with get_hdf5file_replace(outfilename) as o:
    for path in paths:
      try:
        f = get_hdf5file(join(path, infilename))
      except OSError:
        print('No file ', infilename, ' found in ', path)
      else:
        if only_base_path:
          name_ = '/' + basename(normpath(path))
        else:
          name_ = '/' + path

        if source_base_path:
          source_ = '/' + basename(normpath(path))
        else:
          source_ = '/'

        f.copy(source=source_, dest=o, name=name_)
        f.close()

def copy_hdf5_from_subfolders_to_single_file(path, infilename: str, outfilename: str, only_base_path: bool = False, source_base_path: bool = False):
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
  only_base_path, source_base_path: passed to
    copy_from_paths_to_single_file, please look there for further
    information.
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

  copy_hdf5_from_paths_to_single_file(folders, infilename, outfilename, only_base_path, source_base_path)


def clean_up_after_run(tag_first:int, tag_last: int):
  import os

  # Delete single HDF5 files
  os.remove('propagator_0_0.h5')
  os.remove('taginfo.h5')

  # Delete indexed HDF5 files
  for k in range(tag_first, tag_last+1):
    integer_part = '{0}'.format(k).strip()

    try:
      os.remove('spitf_' + integer_part + '.h5')
    except FileNotFoundError:
      pass

    try:
      os.remove('enetf_' + integer_part + '.h5')
    except FileNotFoundError:
      pass

    try:
      os.remove('dentf_' + integer_part + '.h5')
    except FileNotFoundError:
      pass

    try:
      os.remove('phi_mesh_' + integer_part + '.h5')
    except FileNotFoundError:
      pass

    try:
      os.remove('sizeplot_etalev_' + integer_part + '.h5')
    except FileNotFoundError:
      pass

    for l in range(tag_first, tag_last+1):
      two_integer_part = '{0}_{1}'.format(k,l).strip()
      try:
        os.remove('propagator_' + two_integer_part + '.h5')
      except FileNotFoundError:
        pass

      try:
        os.remove('propagator_boundary_' + two_integer_part + '.h5')
      except FileNotFoundError:
        pass

      try:
        os.remove('reconstruct_' + two_integer_part + '.h5')
      except FileNotFoundError:
        pass

      try:
        os.remove('binarysplit_' + two_integer_part + '.h5')
      except FileNotFoundError:
        pass


def prop_reconstruct_3(outfilename: str= 'final.h5'):
  """Collect results for PAR version of neo-2 as function of same name.

  This function is intended to copy the functionality of the subroutine
  prop_reconstruct_3 of the PAR version of neo-2.
  This means it will collect the (most) hdf5 output in a single hdf5
  file, and delete most of the collected files (efinal.h5, fulltransp.h5
  and neo2_config.h5 are the exceptions).

  So far, compared to this subroutine there is an addition and a
  difference.

  Additionaly the name of the output file can be selected.
  The difference is, that this function does not use the 'lsw_save_'
  switches. It simply checks if the files are present, and if so adds
  the content to the output file and deletes them afterwards.

  input:
  ------
  outfilename: str, name(+path) of the output file. Defaults to
    'final.h5'.
  """
  import os

  print("Merging HDF5 files...\n")

  #**********************************************************
  # Open taginfo
  #**********************************************************
  with get_hdf5file('taginfo.h5') as taginfo:
    tag_first = taginfo['tag_first'].value
    tag_last = taginfo['tag_last'].value

  with get_hdf5file('neo2_config.h5') as config:
    if (config['/collision/isw_axisymm'].value == 1):
      tag_first = 3
      tag_last = 3

  #**********************************************************
  # Get current working directory
  #**********************************************************
  cwd = os.getcwd()
  surfname = os.path.basename(cwd).strip()
  print("Using " + surfname + " as surfname.")

  #**********************************************************
  # Create result file
  #**********************************************************
  final =  get_hdf5file_new(outfilename)
  surf = final.create_group(surfname)
  neo2 = surf.create_group('NEO-2')
  propagators = neo2.create_group('propagators')

  cg0_1_num_prop = []
  cg2_1_num_prop = []
  cg0_2_num_prop = []
  cg2_2_num_prop = []
  cg0_3_num_prop = []
  cg2_3_num_prop = []
  denom_mflint_prop = []

  #**********************************************************
  # Iterate over all propagators and merge files
  #**********************************************************
  for k in range(tag_first, tag_last+1):
    h5_filename = '{0}'.format(k).strip()
    prop = propagators.create_group(h5_filename)

    if os.path.exists('spitf_' + h5_filename + '.h5'):
      with get_hdf5file('spitf_' + h5_filename + '.h5') as propfile:
        propfile.copy(source='/', dest=prop, name='spitf')

    if os.path.exists('dentf_' + h5_filename + '.h5'):
      with get_hdf5file('dentf_' + h5_filename + '.h5') as propfile:
        propfile.copy(source='/', dest=prop, name='dentf')

    if os.path.exists('enetf_' + h5_filename + '.h5'):
      with get_hdf5file('enetf_' + h5_filename + '.h5') as propfile:
        propfile.copy(source='/', dest=prop, name='enetf')

    with get_hdf5file('phi_mesh_' + h5_filename + '.h5') as propfile:
      cg0_1_num_prop.append(propfile['cg0_1_num'].value)
      cg2_1_num_prop.append(propfile['cg2_1_num'].value)
      cg0_2_num_prop.append(propfile['cg0_2_num'].value)
      cg2_2_num_prop.append(propfile['cg2_2_num'].value)
      cg0_3_num_prop.append(propfile['cg0_3_num'].value)
      cg2_3_num_prop.append(propfile['cg2_1_num'].value)
      denom_mflint_prop.append(propfile['denom_mflint'].value)
      propfile.copy(source='/', dest=prop, name='phi_mesh')

    with get_hdf5file('sizeplot_etalev_'+h5_filename+'.h5') as propfile:
      propfile.copy(source='/', dest=prop, name='sizeplot_etalev')

  cg0_1_avg = sum(cg0_1_num_prop) / sum(denom_mflint_prop)
  cg2_1_avg = sum(cg2_1_num_prop) / sum(denom_mflint_prop)
  cg0_2_avg = sum(cg0_2_num_prop) / sum(denom_mflint_prop)
  cg2_2_avg = sum(cg2_2_num_prop) / sum(denom_mflint_prop)
  cg0_3_avg = sum(cg0_3_num_prop) / sum(denom_mflint_prop)
  cg2_3_avg = sum(cg2_3_num_prop) / sum(denom_mflint_prop)

  print("cg0_1 = {}".format(cg0_1_avg))
  print("cg2_1 = {}".format(cg2_1_avg))
  print("cg0_2 = {}".format(cg0_2_avg))
  print("cg2_2 = {}".format(cg2_2_avg))
  print("cg0_3 = {}".format(cg0_3_avg))
  print("cg2_3 = {}".format(cg2_3_avg))

  #**********************************************************
  # Merge additional files
  #**********************************************************
  with get_hdf5file('efinal.h5') as propfile:
    propfile.copy(source='/', dest=neo2, name='efinal')

  with get_hdf5file('fulltransp.h5') as propfile:
    propfile.copy(source='/', dest=neo2, name='fulltransp')

  with get_hdf5file('neo2_config.h5') as propfile:
    propfile.copy(source='/', dest=neo2, name='neo2_config')

  with get_hdf5file('taginfo.h5') as propfile:
    propfile.copy(source='/', dest=neo2, name='taginfo')

    neo2['taginfo/cg0_1_avg'] = cg0_1_avg
    neo2['taginfo/cg2_1_avg'] = cg2_1_avg

    neo2['taginfo/cg0_2_avg'] = cg0_2_avg
    neo2['taginfo/cg2_2_avg'] = cg2_2_avg

    neo2['taginfo/cg0_3_avg'] = cg0_3_avg
    neo2['taginfo/cg2_3_avg'] = cg2_3_avg

  print('Result file created. Deleting old files....')

  clean_up_after_run(tag_first, tag_last)

  print("Done.\n")


def merge_maximum(a: list, b: list):
  """Return a merged list by taking the maximum of each entry.

  Example:
  l = merge_maximum([2, 3, 4, 5], [3, 2, 1, 6])

  would mean that
  l = [3, 3, 4, 6]

  input:
  ------
  a, b: lists, which to merge. Should have the same number of elements,
    but this is not checked by this function.
  """
  return [max(*z) for z in zip(a,b)]


def create_higher_dimensional_groupstructure(source, destination, size: int, ignore_version: bool = True):
  for k in source.keys():
    if (k.lower() == 'version' and ignore_version):
      continue
    # ~ if (type(source[k]) == type(h5py.group)):
      # ~ sub = destination.create_group(k)
      # ~ create_higher_dimensional_groupstructure(k, sub, size)
    # ~ else:
    shape = list(source[k].shape)
    shape.insert(0, size)
    destination.create_dataset(k, tuple(shape), dtype=source[k].dtype)


def get_types_of_keys_in_group(group, types: dict, ignore_version: bool = True):
  """Return the types of the keys in a given group as dictionary.

  The function is used recursively for subgroups.

  input:
  ------
  group: h5py group, for which to create the dictionary.
  types: dict, dictionary to which to add the keys. For parent group of
    hdf5 file you probably want to pass an empty dictionary.
  ignore_version: bool, if True fields named 'version' are ignored for
    creating the dictionary. [True]

  output:
  -------
  (recursive) dictionare where each entry is the type of the
  corresponding key.
  """
  for k in group.keys():
    if type(group[k]) == type(group):
      types[k] = {}
      types[k] = get_types_of_keys_in_group(group[k], types[k], ignore_version)
    else:
      types[k] = group[k].dtype

  return types


def get_types_of_keys(name_part: str, tag_first: int, ignore_version: bool = True):
  """Return the types that exist in a file as dictionary.

  Assumes specific format of filename:

  name_part + _ + #tag_first + '.h5'

  Where #tag_first means the string representation of the input
  parameter.

  input:
  ------
  name_part: str, name part of the file from which to extract types.
  tag_first: int, number part of the file from which to extract types.
  ignore_version: bool, if True fields named 'version' are ignored for
    creating the dictionary. [True]

  output:
  -------
  (recursive) dictionare where each entry is the type of the
  corresponding key.
  """
  types = {}
  with get_hdf5file(name_part + '_{0}'.format(tag_first).strip() + '.h5') as f:
    types = get_types_of_keys_in_group(f, types, ignore_version)

  return types


def get_size_of_keys_in_group(f, d: dict, ignore_version: bool = True):
  """Determine size of data elements in group and return it as dictionary.

  This function is used recursively for subgroups.

  input:
  ------
  f: h5py group, for which to determine sizes.
  d: dictionary, if no entry for a dataset exists, it is created. If one
    exists already, then for each dimension the maximum is used. If this
    function is called for the root '/', then this should be an empty
    dictionary, unless you want to make sure that fields have a minimum
    size.
  ignore_version: bool, if True fields named 'version' are ignored for
    creating the dictionary. [True]
  """
  for k in f.keys():
    if (k.lower() == 'version' and ignore_version):
      continue
    if type(f[k]) == type(f):
      if k not in d:
        d[k] = {}
      d[k] = ghi(f[k], d[k], ignore_version)
    else:
      if k not in d:
        d[k] = list(f[k].shape)
      else:
        d[k] = merge_maximum(d[k], list(f[k].shape))

  return d


def get_maximum_sizes_of_groups(name_part: str, tag_first: int, tag_last: int, ignore_version:bool = True):
  """Get maximum size of the groups in a batch of files as dictionary.

  Get the maximum size of datasets within a batch of files, that follow
  the naming convention

  name_part_#.h5

  where '#' stands for the string representation of a number from
  tag_first to tag_last (including).

  input:
  ------
  name_part: str, name part of the file from which to extract types.
  tag_first, tag_last: integers, giving the range of the number part of
    the file from which to extract sizes.
  ignore_version: bool, if True fields named 'version' are ignored for
    creating the dictionary. [True]

  output:
  -------
  Dictionary, where entry of a key is either a list which corresponds to
  the maximum size of that entry in the files, or a dictionary if the
  key corresponds to a group.
  """
  d = {}
  for t in range(tag_first, tag_last+1):
    integer_part = '_{0}'.format(t).strip()
    with get_hdf5file(name_part + integer_part + '.h5') as f:
      get_size_of_keys_in_group(f, d, ignore_version)

  return d


def construct_group(group, sizes: dict, types: dict, variable_length_array:bool = False):
  """ Create subgroups for a given h5py group.

  This function will create a hierachie of subgroups recursively, by
  using sizes and types given as dictionaries.

  input:
  ------
  group: h5py group
  sizes: dictionary, if a key corresponds to a dataset it contains a
    list of integers giving the sizes for each dimension. If it
    corresponds to a group, then it contains a dictionary with the same
    structure.
  types: dictionary, same structure as sizes, but having h5py dtypes as
    entries.
  variable_length_array: bool, if true then two dimensional datasets are
    created as variable length arrays, i.e. only the first dimension is
    fixed. [False]
  """
  import h5py

  for k in sizes.keys():
    if type(sizes[k]) == type(dict()):
      subgroup = group.create_group(k)
      construct_group(subgroup, sizes[k], types[k])
    else:
      dt = types[k]
      size = tuple(sizes[k])
      if variable_length_array and len(size) == 2:
        # required for later versions (1.12+)?
        # ~ dt = h5py.vlen_dtype(types[k])
        dt = h5py.special_dtype(vlen=types[k])
        size = (size[0],)
      group.create_dataset(k, size, dt)


def insert_size(sizes: dict, size_to_insert: int):
  """Prepend size to entries in a dictionary.

  This works on dictionaries thise created by get_maximum_sizes_of_groups
  and prepends a size for all entries recursively.


  input:
  ------
  sizes: dictionary, in all entries the given size is inserted as first
    entry.
  size_to_insert: integer, size which to insert.
  """
  for k in sizes.keys():
    if type(sizes[k]) == type(dict()):
      insert_size(sizes[k], size_to_insert)
    else:
      sizes[k].insert(0, size_to_insert)

  return sizes


def copy_group(source, destination, index: int, ignore_version:bool = True, variable_length_array:bool = False):
  """ Copy the content of one group to another recursively.

  Copy the content of one file/group to another, with two assumptions
  - destination is already created with additional first dimension.
  - copying is done to the other dimension, which are at least large
    enough to hold the elements, but might be bigger.

  input:
  ------
  source: h5py group from which to copy data.
  destination: h5py group to which to copy data.
  index: integer, first index for entries in destination.
  ignore_version: bool, if True fields named 'version' are ignored for
    creating the dictionary. [True]
  variable_length_array: bool, if this is true then one dimensional
    arrays from the source, are stored as variable length arrays in
    destination. [False]
  """
  for k in source.keys():
    if (k.lower() == 'version' and ignore_version):
      continue
    if type(source[k]) == type(dict()):
      copy_group(source[k], destionation[k], index, ignore_version)
    else:
      size = source[k].shape
      try:
        if len(size) == 0:
          destination[k][index] = source[k].value
        elif len(size) == 1:
          if variable_length_array:
            destination[k][index] = source[k]
          else:
            destination[k][index, 0:size[0]] = source[k]
        elif len(size) == 2:
          destination[k][index, 0:size[0], 0:size[1]] = source[k]
        elif len(size) == 3:
          destination[k][index, 0:size[0], 0:size[1], 0:size[2]] = source[k]
        elif len(size) == 4:
          destination[k][index, 0:size[0], 0:size[1], 0:size[2], 0:size[3]] = source[k]
        elif len(size) == 5:
          destination[k][index, 0:size[0], 0:size[1], 0:size[2], 0:size[3], 0:size[4]] = source[k]
        elif len(size) == 6:
          destination[k][index, 0:size[0], 0:size[1], 0:size[2], 0:size[3], 0:size[4], 0:size[5]] = source[k]
        else:
          raise Exception('copy_group: dimension of source array not yet implemented!')
      except KeyError:
        print('Key error in copy_group: len(sizes) = {0}, key = {1}, index = {2}'.format(len(sizes), k, index))
        raise
      except TypeError:
        print('Type error in copy_group: sizes = {0}, shape of source = {1}'.format(list(sizes[k]), source[k].shape))
        raise


def prop_reconstruct_3a(outfilename: str= 'testing_only_final.h5', variable_length_array:bool = False):
  """
  input:
  ------
  outfilename: str, name(+path) of the output file. Defaults to
    'final.h5'.
  variable_length_array: bool, if this is true, then one dimensional
   arrays from the propagator inputs will be stored as variable length
   arrays.
  """
  import os

  ignore_version_ = True

  print("Merging HDF5 files...\n")

  #**********************************************************
  # Open taginfo
  #**********************************************************
  with get_hdf5file('taginfo.h5') as taginfo:
    tag_first = taginfo['tag_first'].value
    tag_last = taginfo['tag_last'].value

  with get_hdf5file('neo2_config.h5') as config:
    if (config['/collision/isw_axisymm'].value == 1):
      tag_first = 3
      tag_last = 3

  #**********************************************************
  # Get current working directory
  #**********************************************************
  cwd = os.getcwd()
  surfname = os.path.basename(cwd).strip()
  print("Using " + surfname + " as surfname.")

  #**********************************************************
  # Create result file
  #**********************************************************
  final = get_hdf5file_new(outfilename)
  surf = final.create_group(surfname)

  cg0_1_num_prop = []
  cg2_1_num_prop = []
  cg0_2_num_prop = []
  cg2_2_num_prop = []
  cg0_3_num_prop = []
  cg2_3_num_prop = []
  denom_mflint_prop = []

  file_batches_to_consider = ['spitf', 'dentf', 'enetf', 'phi_mesh', 'sizeplot_etalev']
  for f in file_batches_to_consider:
    if os.path.exists(f + '_{0}'.format(tag_first) + '.h5'):
      types = get_types_of_keys(f, tag_first, ignore_version=ignore_version_)
      sizes = get_maximum_sizes_of_groups(f, tag_first, tag_last, ignore_version=ignore_version_)

      #**********************************************************
      # Create data structure for result file
      #**********************************************************
      insert_size(sizes, size_to_insert = tag_last - tag_first + 1)
      destination = surf.create_group(f)
      construct_group(destination, sizes, types, variable_length_array=variable_length_array)

      #**********************************************************
      # Copy the data
      #**********************************************************
      for t in range(tag_first, tag_last+1):
        integer_part = '_{0}'.format(t).strip()

        with get_hdf5file(f + integer_part + '.h5') as source:
          copy_group(source, destination, index= t - tag_first, ignore_version=ignore_version_, variable_length_array=variable_length_array)

          if f == 'phi_mesh': # ~ with get_hdf5file('phi_mesh_' + h5_filename + '.h5') as propfile:
            cg0_1_num_prop.append(source['cg0_1_num'].value)
            cg2_1_num_prop.append(source['cg2_1_num'].value)
            cg0_2_num_prop.append(source['cg0_2_num'].value)
            cg2_2_num_prop.append(source['cg2_2_num'].value)
            cg0_3_num_prop.append(source['cg0_3_num'].value)
            cg2_3_num_prop.append(source['cg2_1_num'].value)
            denom_mflint_prop.append(source['denom_mflint'].value)

  cg0_1_avg = sum(cg0_1_num_prop) / sum(denom_mflint_prop)
  cg2_1_avg = sum(cg2_1_num_prop) / sum(denom_mflint_prop)
  cg0_2_avg = sum(cg0_2_num_prop) / sum(denom_mflint_prop)
  cg2_2_avg = sum(cg2_2_num_prop) / sum(denom_mflint_prop)
  cg0_3_avg = sum(cg0_3_num_prop) / sum(denom_mflint_prop)
  cg2_3_avg = sum(cg2_3_num_prop) / sum(denom_mflint_prop)

  print("cg0_1 = {}".format(cg0_1_avg))
  print("cg2_1 = {}".format(cg2_1_avg))
  print("cg0_2 = {}".format(cg0_2_avg))
  print("cg2_2 = {}".format(cg2_2_avg))
  print("cg0_3 = {}".format(cg0_3_avg))
  print("cg2_3 = {}".format(cg2_3_avg))

  #**********************************************************
  # Merge additional files
  #**********************************************************
  with get_hdf5file('efinal.h5') as propfile:
    propfile.copy(source='/', dest=surf, name='efinal')

  with get_hdf5file('fulltransp.h5') as propfile:
    propfile.copy(source='/', dest=surf, name='fulltransp')

  with get_hdf5file('neo2_config.h5') as propfile:
    propfile.copy(source='/', dest=surf, name='neo2_config')

  with get_hdf5file('taginfo.h5') as propfile:
    propfile.copy(source='/', dest=surf, name='taginfo')

    surf['taginfo/cg0_1_avg'] = cg0_1_avg
    surf['taginfo/cg2_1_avg'] = cg2_1_avg

    surf['taginfo/cg0_2_avg'] = cg0_2_avg
    surf['taginfo/cg2_2_avg'] = cg2_2_avg

    surf['taginfo/cg0_3_avg'] = cg0_3_avg
    surf['taginfo/cg2_3_avg'] = cg2_3_avg

  print('Result file created. Deleting old files....')

  clean_up_after_run(tag_start, tag_last)

  print("Done.\n")


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
  try:
    h5o = get_hdf5file(other_filename)

  except OSError as e:
    print('ERROR: hdf5-Data file "' + other_filename + '" does not exist in test-folder!')
    return [False, False]

  else:
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
