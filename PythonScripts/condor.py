#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 01 13:18:41 2019

@author: Rico Buchholz
"""

class CondorSubmitFile:
  """Class representing a file used by condor when submitting jobs.

  This class represents the information in a file used to commit jobs
  via condor.
  It will parse an existing file to fill the parameters and
  can be written to a file.

  You can get a list of parameters, and the value of a parameter.

  There are methods for setting some of the parameters.

  The two functions condor_submit_file_to_dict and
  dict_to_condor_submit_file do not require an object.
  """


  def condor_submit_file_to_dict(filename: str):
    """Create a dictionary from a htcondor submit file.

    Example usage:
    [d, q] = condor_submit_file_to_dict('./submit')

    d will be the dictionary with parameters and q the line(s) containing
    the queue command.
    As queue lines are considered lines (either of these):
    - which first element is 'queue'
    - which contain only one element
    - whose second element is not '='

    input:
    ------
    filename: str, name(+path) of the file to load.

    output:
    -------
    list with two elements. First is a dictionary with the parameter
    entries of the file. The keys are converted to lowercase.
    The second entry is a string holding the queue command line, (which is
    not included in the dict).
    """
    with open(filename) as f:
      lines = f.readlines()

    lines = [l.strip() for l in lines]

    d = dict()
    q = ''
    for l in lines:
      if len(l) > 0:
        if l[0] == '#':
          continue
        elif l.split()[0].lower() == 'queue':
          q += l
        elif len(l.split()) == 1:
          q += l
        elif l.split()[1] != '=':
          q += l
        else:
          d[l.split('=')[0].strip().lower()] = ''.join(l.strip('\n').split('=', maxsplit=1)[1].strip())

    q = q.replace('\\', ' ')

    return [d, q]


  def dict_to_condor_submit_file(filename: str, d: dict, q: str, overwrite:bool = False):
    """Create a condor submit file from dictionary and queue string.

    Counterpart to condor_submit_file_to_dict.
    Takes a dictionary and a queue string and writes them to a text file,
    that can be used with 'condor_submit'.

    input:
    ------
    filename: string, name (and path) of the file which to write.
    d: dict, parameter values pairs that should be written to file.
    q: string, queue command to be written to file.
    overwrite: bool, if false, an exception will be raised if 'filename'
      already exists. If true, existing files will be overwritten.

    output:
    -------
    None.

    side effects:
    -------------
    Creates a file. If 'overwrite' is true, an existing file will be
    replaced.
    """
    from os import path
    if (path.exists(filename) and not overwrite):
      raise FileExistsError('File already exist and should not be overwriten: ' + filename)

    with open(filename, 'w') as f:
      for k in d.keys():
        f.write(k + ' = ' + d[k] + '\n')

      f.write('\n')

      f.write(q + '\n')


  def read(self, filename: str):
    """Reads in the file.

    This method will open the file and read its content. This is then
    parsed for keywords to set the member of this class.

    input
    -----
    filename: name of the file which should be read.

    TODO
    - Handle the case of non-existing file.
    """

    [self._parameters, self._queue_command] = CondorSubmitFile.condor_submit_file_to_dict(filename)


  def write(self, filename: str, overwrite: bool = False):
    """Create a suitable condor file from this object.

    input
    -----
    filename: name under which to save the file.
    """
    CondorSubmitFile.dict_to_condor_submit_file(filename, self._parameters, self._queue_command, overwrite)


  def set_memory(self, value: int):
    """ Set memory requested by the job(s).

    input:
    ------
    value: integer, memory that should be requested in GB.
    """
    self._parameters['request_memory'] = '{0}*1024'.format(value)


  def _set_cpus_in_arguments(self, value: int):
    """Set number of cpus in argument parameter.

    This functions assumes the 'argument' parameter exists already and
    has a key 'arguments', which contains a part '-np #', where here
    '#' stands for a positive interger.

    If the 'argument' parameter does not exit, then nothing is done.

    input:
    ------
    value: integer, the number of cpus to request.
    """
    if 'arguments' in self._parameters.keys():
      elements = self._parameters['arguments'].split()
      if '-np' in elements:
        for i,x in enumerate(elements):
          if x == '-np':
            elements[i+1] = '{0}'.format(value)
            break
        self._parameters['arguments'] = ' '.join(elements)


  def set_cpus(self, value: int):
    """
    input:
    ------
    value: integer, the number of cpus to request.
    """
    self._parameters['request_cpus'] = '{0}'.format(value)
    self._set_cpus_in_arguments(value)


  def set_requirements(self, value: str):
    """ Set requirements.

    \note: this function does not check if the given string is a valid
      list of requirements.
    """
    self._parameters['requirements'] = value


  def add_requirements(self, value: str):
    """ Add requirement with and.

    \note: this function does not check if the given string is a valid
      requirement.

    input:
    ------
    value: str, containing the requirement, e.g.
      'OpSysAndVer == "Debian10"' or
      'TARGET.UtsnameNodename != "faepop27"'.
    """
    if 'requirements' in self._parameters.keys():
      self._parameters['requirements'] = ' '.join(self._parameters['requirements'].split()[0:-1]) + ' && ' + value + ' )'
    else:
      self._parameters['requirements'] = '( ' + value + ' )'


  def get_keys(self):
    """Return list of parameters known.

    Does not include queue command.
    """
    return list(self._parameters.keys())


  def get_value(self, key: str):
    """Return value of given parameter.

    Does not include queue command.

    input:
    ------
    key: string, name of the parameter for which to return the value.
    """
    return self._parameters[key]


  def __init__(self, filename: str):
    self.read(filename)
