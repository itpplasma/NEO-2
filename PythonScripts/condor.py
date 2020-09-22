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
  It will parse an existing file to fill the parameters (it knows) and
  can be written to a file.

  There are methods for setting some of the parameters.

  TODO:
  The queue_command is so far not treated correctly. It needs to be made
  a distinction between using name expansion (e.g. 'es_*') and using a
  list of files. This is not implemented yet.

  Maybe use htcondor.Submit instead?
  """

  def read(self, filename: str):
    """Reads in the file.

    This method will open the file and read its content. This is then
    parsed for keywords to set the member of this class.

    input
    -----
    filename: name of the file which should be read.

    TODO
    - Handle the case of non-existing file.
    - The queue_command and the folder paramerts are not set correctly
      so far.
    """
    with open(filename) as f:
      lines = f.readlines()

    _Executable = '/usr/bin/mpiexec'
    _arguments = '-mca orte_tmpdir_base \"/tmp/\" -np  2 ./neo_2.x'
    _Universe = 'vanilla'
    _Error      = '../err.$(Process)'
    _Log        = '../scan.log'
    _run_as_owner = True
    _notify_user  = 'buchholz@tugraz.at'
    _notification = 'Error'
    _nice_user    = False
    _requirements = ''
    _request_cpus = 1
    _request_memory = 21
    _Getenv = True

    _Output     = 'out'
    _initialdir = '$(dirname)'
    _queue_command = 'Queue dirname matching dirs'
    _folders = 'es_*'

    for line in lines:
      if (len(line.strip()) > 0):
        if (line.split('=', 1)[0].strip() == 'Executable'):
          self._Executable = line.split('=', 1)[1].strip()
        if (line.split('=', 1)[0].strip() == 'arguments'):
          self._arguments = line.split('=', 1)[1].strip()
        if (line.split('=', 1)[0].strip() == 'Universe'):
          self._Universe = line.split('=', 1)[1].strip()
        if (line.split('=', 1)[0].strip() == 'Error'):
          self._Error = line.split('=', 1)[1].strip()
        if (line.split('=', 1)[0].strip() == 'Log'):
          self._Log = line.split('=', 1)[1].strip()
        if (line.split('=', 1)[0].strip() == 'run_as_owner'):
          self._run_as_owner = bool(line.split('=', 1)[1].strip())
        if (line.split('=', 1)[0].strip() == 'notify_user'):
          self._notify_user = line.split('=', 1)[1].strip()
        if (line.split('=', 1)[0].strip() == 'notification'):
          self._notification = line.split('=', 1)[1].strip()
        if (line.split('=', 1)[0].strip() == 'nice_user'):
          self._nice_user = bool(line.split('=', 1)[1].strip())
        if (line.split('=', 1)[0].strip() == 'requirements'):
          self._requirements = line.split('=', 1)[1].strip()
        if (line.split('=', 1)[0].strip() == 'request_cpus'):
          self._request_cpus = int(line.split('=', 1)[1].strip())
        if (line.split('=', 1)[0].strip() == 'request_memory'):
          self._request_memory = int(line.split('=', 1)[1].strip().split('*')[0])
        if (line.split('=', 1)[0].strip() == 'GetEnv'):
          self._Getenv = bool(line.split('=', 1)[1].strip())
        if (line.split('=', 1)[0].strip() == 'Output'):
          self._Output = line.split('=', 1)[1].strip()
        if (line.split('=', 1)[0].strip() == 'initaldir'):
          self._initialdir = line.split('=', 1)[1].strip()
        if (line.split()[0].strip() == 'Queue'):
          self._queue_command = line

  def _get_boolean_string(self, value: bool):
    """Simple function to convert a python bool to a condor bool (as string).
    """
    if value:
      boolean_string = 'true'
    else:
      boolean_string = 'false'
    return boolean_string

  def write(self, filename: str):
    """Create a suitable condor file from this object.

    input
    -----
    filename: name under which to save the file.

    TODO
    The queue command needs correct treatment.
    """
    lines = 'Executable = ' + self._Executable + '\n'
    lines += 'arguments = ' + self._arguments + '\n'
    lines += 'Universe = ' + self._Universe + '\n'
    lines += 'Error = ' + self._Error + '\n'
    lines += 'Log = ' + self._Log + '\n'
    lines += 'run_as_owner = ' + self._get_boolean_string(self._run_as_owner) + '\n'
    lines += 'notify_user = ' + self._notify_user + '\n'
    lines += 'notification = ' + self._notification + '\n'
    lines += 'nice_user = ' + self._get_boolean_string(self._nice_user) + '\n'
    lines += 'requirements = ' + self._requirements + '\n'
    lines += 'request_cpus = ' + str(self._request_cpus) + '\n'
    lines += 'request_memory = ' + str(self._request_memory) + '*1024\n'
    lines += 'Getenv = ' + self._get_boolean_string(self._Getenv) + '\n'
    lines += '\n'
    lines += 'Output = ' + self._Output + '\n'
    lines += 'initialdir = ' + self._initialdir + '\n'
    lines += self._queue_command #+ ' ' + self._folders + '\n'

    with open(filename, 'w') as f:
      f.writelines(lines);

  def set_folders(self, value: str):
    self._folders = value

  def set_memory(self, value: int):
    self._request_memory = value

  def set_cpus(self, value: int):
    self._request_cpus = value

  def set_requirements(self, value: str):
    self._requirements = value

  def __init__(self, filename: str):
    self.read(filename)
