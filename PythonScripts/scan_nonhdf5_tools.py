#!/usr/bin/env python3

class condor_run:
  """Class for storing information about one condor job.

  This class is intended for storing information about one condor job.
  This includes scan and job id, submit date and submit host, the host
  on which the job did run, the start date and end date, and the maximum
  memory consumption.

  Dates are stored as datetime object. Note that the condor log does
  not include a year, thus old jobs or jobs run around new year might
  cause problems.
  """

  def parse_id_string(self, id_string: str):
    """Parse an id string from a condor log and return the information.

    Example of id_string: '(1807.000.000)'
    The first part is called here the scan id, the second the run id.
    The third part is ignored so far.
    """
    scan_id = int(id_string[1:-2].split('.')[0])
    run_id = int(id_string[1:-2].split('.')[1])
    return [scan_id, run_id]

  def parse_datetime_string(self, datetime_string: str):
    """Parse a datetime string from a condor log and return the information.

    Example of datetime_string: '04/25 13:28:02'
    Order is month, day, hour, minutes second. Note that the year is not
    provided.
    """
    from datetime import datetime

    return datetime.strptime(datetime_string, '%m/%d %H:%M:%S')

  def __init__(self, text_submit_messgage):
    """Initalize the object.

    Example of message: '000 (1807.000.000) 04/25 13:28:02 Job submitted from host: <129.27.161.38:9618?addrs=129.27.161.38-9618&noUDP&sock=942_b64d_3>'
    The submit message is parsed to get a bunch of information (id's,
    submit date and submit host).
    For the rest, other methods need to be called, they are here only
    set to an empty string.
    """
    parts = text_submit_messgage[0].split()

    [self.scan_id, self.run_id] = self.parse_id_string(parts[1])
    self.submit_date = self.parse_datetime_string(parts[2] + ' ' + parts[3])
    self.submit_host = parts[8]

    self.host = ''
    self.start_date = ''
    self.end_date = ''
    self.num_processors = 1
    self.usertime = 0
    self.systemtime = 0
    self.memory = 0

    self.aborted_by_user = False

  def set_start_date(self, text_start_date: str):
    """Parse passed text and set the start date.
    """
    self.start_date = self.parse_datetime_string(text_start_date)

  def set_end_date(self, text_end_date: str):
    """Parse passed text and set the end date.
    """
    self.end_date = self.parse_datetime_string(text_end_date)

  def set_memory(self, text_memory: str):
    """Parse passed text and set the memory.

    Parsing might be a bit exagerated, the string is just converted to
    an integer.
    """
    self.memory = int(text_memory)

  def set_num_processors(self, text_num_processors: str):
    """Parse passed text and set the number of processors.

    Parsing might be a bit exagerated, the string is just converted to
    an integer.
    """
    self.num_processors = int(text_num_processors)

  def list_time_to_seconds(self, dayshoursminutesseconds):
    """Returns the number of seconds for a list with days:hours:minutes:seconds.
    """
    return dayshoursminutesseconds[0]*24*60*60 + dayshoursminutesseconds[1]*60*60 + dayshoursminutesseconds[2]*60 + dayshoursminutesseconds[3]

  def set_usertime(self, text_usertime: str):
    """Parse passed text and set the usertime.

    Time is expected as a string with days:hours:minutes:seconds.
    """
    if (text_usertime[-1] == ','):
      text_usertime = text_usertime[0:-1]
    try:
      parts = [int(x) for x in text_usertime.split(':')]
    except ValueError:
      print('Excpetion while parsing string for usertime:')
      print(text_usertime)
      raise

    self.usertime = self.list_time_to_seconds(parts)

  def set_systemtime(self, text_systemtime: str):
    """Parse passed text and set the usertime.

    Time is expected as a string with days:hours:minutes:seconds.
    """
    parts = [int(x) for x in text_systemtime.split(':')]
    self.systemtime = self.list_time_to_seconds(parts)

  def parse_start_message(self, text_start_message):
    """Parse the start message to extract information and set members.

    Example of message: 001 (1873.075.000) 06/04 14:15:58 Job executing on host: <129.27.161.108:9618?addrs=129.27.161.108-9618&noUDP&sock=962_6e32_4>
    """
    parts = text_start_message[0].split()
    self.set_start_date(parts[2] + ' ' + parts[3])
    self.host = parts[8]

  def parse_end_message(self, text_end_message):
    """Parse the start message to extract information and set members.

    This includes messages with the identifiers 4 (evicted, e.g. abort
    by user) and 5 (terminated, e.g. finished)

    Example of message:
    005 (1873.074.000) 06/04 14:30:28 Job terminated.
    (1) Normal termination (return value 0)
      Usr 0 00:53:52, Sys 0 00:16:23  -  Run Remote Usage
      Usr 0 00:00:00, Sys 0 00:00:00  -  Run Local Usage
      Usr 0 00:53:52, Sys 0 00:16:23  -  Total Remote Usage
      Usr 0 00:00:00, Sys 0 00:00:00  -  Total Local Usage
    0  -  Run Bytes Sent By Job
    0  -  Run Bytes Received By Job
    0  -  Total Bytes Sent By Job
    0  -  Total Bytes Received By Job
    Partitionable Resources :    Usage  Request Allocated
      Cpus                 :                 2         2
      Disk (KB)            :       75       75    367589
      Memory (MB)          :    25659    31744     31744
    """
    parts = text_end_message[0].split()
    time_index = 2
    if (text_end_message[1].split()[-1] == '7)'):
      # At least on line added. Might be more, then a search has to be
      # done, but this is not implemented yet. (Only user and system
      # time should be affected by this.)
      time_index += 1
    try:
      self.set_end_date(parts[2] + ' ' + parts[3])
      self.set_memory(text_end_message[-1].split()[3])
      self.set_num_processors(text_end_message[-3].split()[3])
      self.set_usertime(text_end_message[time_index].split()[1] + ':' + text_end_message[time_index].split()[2])
      self.set_systemtime(text_end_message[time_index].split()[4] + ':' + text_end_message[time_index].split()[5])
    except ValueError:
      print('Exception while parsing end message:')
      print(text_end_message)
      raise

  def parse_abort_message(self, text_abort_message):
    """Parse abort message.

    This function parses the abort message (identifier 9). So far the
    only action is to set status aborted by user to true.
    \Attention It is not known if this message is also send for other
      reasons, like abort by the system.
    """
    self.aborted_by_user = True

  def get_scan_id(self):
    return self.scan_id

  def get_run_id(self):
    return self.run_id

class condor_log:
  """Class for scanning of the log file of a condor run and storing basic information.

  This is done by creating a list of condor_run objects.

  Example usage:
  # Creating
  l = condor_log('path/filename')

  # Get time needed per job
  t = l.get_time_per_job()

  # Get memory needed per job
  m = l.get_memory_per_job()

  # Plotting, assumes 'plot' is defined.
  # This will plot 'minor' job number on x, and memory on y axis.
  plot([k[1] for k in m[0]], m[1])
  """

  def split_into_text_messages(self, lines):
    """Takes the lines of a file, and splits it into individual messages.

    This takes a list with the lines of a file, and sorts them into
    lists of lists, where each entry of the outer represents a single
    message, and the inner one are the lines of the message.
    Messages are assumed to be separated by lines containing only three
    dots '...'.
    """
    text_messages = []
    text_messages.append([])
    for l in lines:
      if l.strip() != '...':
        text_messages[-1].append(l.strip())
      else:
        text_messages.append([])

    # Messages end with '...' so we added one element to much.
    text_messages.pop()

    return text_messages

  def get_run(self, string_id: str):
    """Get a specific run from the list, with an id given as string.
    """
    [scan_id, run_id] = self.runs[0].parse_id_string(string_id)
    for r in self.runs:
      if (r.get_scan_id() == scan_id) and (r.get_run_id() == run_id):
        return r

  def parse_text_message(self, text_message):
    """Parse a single message.

    This is done by determining the type of the message and call the
    corresponding methods of the condor_run object.
    If a job is submitted, a condor_run is appended.
    If a job is started, the start message is parsed.
    If a job is finished, the end message is parsed.

    The type of the message is determined by the first number in the
    first line:
    0 submit
    1 start
    4 evicted
    5 end
    6 update image size (not used at the moment)
    9 aborted by user
    """
    message_identifier = int(text_message[0].split()[0])

    # submitted
    if message_identifier == 0:
      self.runs.append(condor_run(text_message))
    # job started
    elif message_identifier == 1:
      r = self.get_run(text_message[0].split()[1])
      r.parse_start_message(text_message)
    # job terminated or got evicted
    elif message_identifier == 5 or message_identifier == 4:
      r = self.get_run(text_message[0].split()[1])
      r.parse_end_message(text_message)
    # update of image size.
    elif message_identifier == 6:
      # Do nothing, memory is parsed from end message.
      # Might be interesting in the future to get memory over time.
      None
    elif message_identifier == 9:
      r = self.get_run(text_message[0].split()[1])
      r.parse_abort_message(text_message)

  def parse_text_messages(self, text_messages):
    """Parse the message: go through the list and parse each message.
    """
    for message in text_messages:
      self.parse_text_message(message)

  def __init__(self, filename: str):
    """Opens the file, reads the content, splits and parses it.
    """
    with open(filename) as f:
      self.lines = f.readlines()

    self.runs = []

    self.text_messages = self.split_into_text_messages(self.lines)

    self.parse_text_messages(self.text_messages)

  def get_memory_per_job(self):
    """Return lists with ids and memory consumption (MBy) of the jobs.

    This will return two lists. The first contains the scan and run id
    of the jobs, the second the respective memory consumption of the
    job.
    """
    ids = []
    memory_per_job = []

    for r in self.runs:
      ids.append([r.get_scan_id(), r.get_run_id()])
      memory_per_job.append(r.memory)

    return [ids, memory_per_job]

  def get_max_memory(self):
    """Return maximum memory consumption (MBy).
    """
    [ids, memory_per_job] = self.get_memory_per_job()

    return max(memory_per_job)

  def get_time_per_job(self):
    """Return lists with ids and time (datetime.timedelta) of the jobs.

    This will return two lists. The first contains the scan and run id
    of the jobs, the second the respective time needed of the job, as a
    datetime.timedelta object.
    If an entry could not be set due to a TypeError, then the timedelta
    is set to -1 day.
    """
    from datetime import timedelta

    ids = []
    time_per_job = []

    for r in self.runs:
      ids.append([r.get_scan_id(), r.get_run_id()])
      try:
        time_per_job.append(r.end_date - r.start_date)
      except TypeError:
        # Assumed that this happened because either end_date or end_date
        # and start_date are not set, for example because the job was
        # aborted by the user.
        time_per_job.append(timedelta(days=-1))
      except:
        print('Exception while determining time per job:')
        print(ids[-1])
        print(r.start_date)
        print(r.end_date)
        raise

    return [ids, time_per_job]

  def get_max_time(self):
    """Return maximum time as datetime.timedelta.
    """
    [ids, time_per_job] = self.get_time_per_job()

    return max(time_per_job)

  def get_usertime_per_job(self):
    """Return lists with ids and usertime of the jobs.

    This will return two lists. The first contains the scan and run id
    of the jobs, the second the respective usertime needed of the job
    in seconds.
    """
    ids = []
    time_per_job = []

    for r in self.runs:
      ids.append([r.get_scan_id(), r.get_run_id()])
      time_per_job.append(r.usertime)

    return [ids, time_per_job]

  def get_systemtime_per_job(self):
    """Return lists with ids and systemtime of the jobs.

    This will return two lists. The first contains the scan and run id
    of the jobs, the second the respective systemtime needed of the job
    in seconds.
    """
    ids = []
    time_per_job = []

    for r in self.runs:
      ids.append([r.get_scan_id(), r.get_run_id()])
      time_per_job.append(r.systemtime)

    return [ids, time_per_job]


def get_memory_consumption_of_run(lines, id_of_run: str):
  """Get memory consumption for 'id_of_run' from 'lines'.

  Look for the memory consumption of a certain run in a textfile, given
  as list of lines.
  The output for the memory is expected to be in a certain format,
  example:

004 (1708.099.000) 12/13 08:47:08 Job was evicted.
  (0) Job was not checkpointed.
    Usr 0 14:48:28, Sys 0 00:16:42  -  Run Remote Usage
    Usr 0 00:00:00, Sys 0 00:00:00  -  Run Local Usage
  0  -  Run Bytes Sent By Job
  0  -  Run Bytes Received By Job
  Partitionable Resources :    Usage  Request Allocated
     Cpus                 :                 2         2
     Disk (KB)            :       75       75   1627536
     Memory (MB)          :    17797    35840     35840

  Important here are the offset of the line with the memory (the last
  one in this example) from the line with the id of the run (first line
  in this example), and the whitespace in the line with the memory. For
  the whitespace there should be no problem as long as there is at least
  one 'space' between the elements.

  Parameter
  ----------
  lines: list of string, that represent a textfile.
  id_of_run: a string which identifies the run for which to return
    memory consumption. Expected to be in the format number.###, where
    number identifies a scan.

  Result
  ----------
  Integer, giving the memory consumption in megabytes. If no indication
    of memory consumption is found in the lines, then -1 is returned.
  """

  run_finished_phrase='Job terminated.'
  offset_line_with_memory_consumption=13
  # Remember
  nr_memory_element_in_line=4
  offset_memory_element_in_line=nr_memory_element_in_line-1
  for i in range(len(lines)):
    if (lines[i].find(id_of_run) >= 0 and lines[i].find(run_finished_phrase) >= 0):
      parts = lines[i+offset_line_with_memory_consumption].split()
      return int(parts[offset_memory_element_in_line])

  return int(-1)

def get_memory_consumption_of_scan(lines, id_scan: str, num_runs: int):
  """Get memory consumption from 'lines' for a scan, with total of 'num_runs'

  This will return a list with the memory consumption of runs belonging
  to a scan.
  This assumes that id_run= id_scan.###.
  It is not sure if this function will work if a scan contained more
  than 1000 runs.

  Makes use of get_memory_consumption_of_run.

  Parameter
  ----------
  lines: list of string, that represent a textfile.
  id_scan: string that contains a number, e.g. '1708', which identifies
    the scan.

  Result
  ----------
  List of integers with 'num_runs' elements, giving the memory
    consumption for each of the runs.
  """

  memory_consumption = []
  for i in range(num_runs):
    memory_consumption.append(get_memory_consumption_of_run(lines, id_scan+".{:03d}".format(i)))

  return memory_consumption

def set_neo2in(folder: str, subfolder_pattern: str, vphifilename: str, backup: bool, read_precom: bool):
  """Set read/write switches and vphi in neo2in files in subfolders.

  This function is intended to set the switches lsw_read_precom and
  lsw_write_precom for a computation run, i.e. the former is set to true
  and the latter to false.
  Also, vphi will be adjusted, the shift will be read from the file
  vphifilename, and if this is not found, a default is used. The file
  must contain the value in the form 'vphiref = 0' in the first line.

  Example usage:
  from scan_nonhdf5_tools import set_neo2in
  set_neo2in('./', 'es_*', 'vphiref.in', True, True)

  input:
  ------
  folder: string with the folder where the subfolders are located. If
    you are in the same directory, then './' can be used.
  subfolder_pattern: string which describes the subfolders. This may
    contain wildcards, e.g. 'es_*'.
    Attention: the pattern should not include files that are also
    present in the directory, as there is no check, to operate on
    folders only.
  vphifilename: string, with the name of the file, that contains the
    value for the shift in vphi. If the file can not be found, this is
    not considered an error, instead a default value is used.
  backup: boolean, if true, a backup of the original input file is
    stored with added '~' at the end.
  """
  from os.path import join
  from pathlib import Path

  import f90nml.parser

  # Check if vphiref.in is present.
  #   if so, extract the value of vphi.
  # create list of folders/iterate over folders.
  # for each folder
  #   read in neo2.in
  #   eventually backup neo2.in
  #   set lsw_read_precom to true and lsw_write_precom to false.
  #   set vphi
  #   write neo2.in

  vphi = 0.0
  try:
    with open(join(folder, vphifilename)) as f:
      line = f.readline()
      parts = line.split('=')
      if (len(parts) > 1):
        vphi = float(parts[1].split('!')[0])
  except FileNotFoundError:
    # If the file is not found, this is not considered to be an error,
    # just leave vphi on the default value.
    print('\n')
    print('File "' + vphifilename + '" not found.')
    print('Using default value: {:e}'.format(vphi))
    print('\n')

  # Get all objects in the folder
  p = Path(folder)
  #~ # Filter for the directories
  #~ p = [x for x in p.iterdir() if x.is_dir()]
  # Filter for the name of the subdirectories
  folders = list(p.glob(subfolder_pattern))

  parser = f90nml.parser.Parser()
  for d in folders:
    current_filename = join(folder, d.name, 'neo2.in')
    nml = parser.read(current_filename)

    print('Processing: ' + current_filename)

    if (backup):
      nml.write(current_filename + '~')

    # Adjust the rotation, new value is background + shift, assuming
    # The file contains already the background value.
    name_value_tuples = [('collision', 'lsw_read_precom', read_precom),
                         ('collision', 'lsw_write_precom', False),
                         ('multi_spec', 'vphi', nml['multi_spec']['vphi'] + vphi)]
    change_namelist_values_for_file_object(nml, name_value_tuples)
    nml.write(current_filename, True)

def set_neo2in_reconstruction(folder: str, subfolder_pattern: str, backup: bool, value: int):
  """
  \brief Top-level function to reset input variable 'prop_reconstruct'.

  This is a top level function intended for use with the stellerator
  variant of neo2.
  It will set the input variable 'prop_reconstruct' to given value for
  all subfolders that match a given pattern in a given folder.

  WARNING: This function does not check if the given value is valid.
    This is intentionaly to avoid having to change this function if the
    range of valid values changes.

  input:
  ------
  folder: string, including location, determining where to find the
    subfolders. If they are in the local folder './' is an appropriate
    value.
  subfolder_pattern: string, determining names of the subfolders where
    to look for neo2.in files to change. Example would be 'es_*' and
    's*'.
  backup: boolean, if true the original version is saved with '~'
    appended to the name.
  value: integer, the value to use for 'prop_reconstruct'.

  output:
  -------
  No formal output.

  Side-effects:
  -------------
  Changes the files determined by input parameters, and may create
  copies of the older version.
  """
  from pathlib import Path

  # Get all objects in the folder
  p = Path(folder)
  #~ # Filter for the directories
  #~ p = [x for x in p.iterdir() if x.is_dir()]
  # Filter for the name of the subdirectories
  folders = list(p.glob(subfolder_pattern))

  set_neo2in_reconstruction_for_folders(folders, folder, backup, value)


def set_neo2in_reconstruction_for_folders(folders, folder: str, backup: bool, value: int):
  """
  \brief Medium-level function to reset input variable 'prop_reconstruct'.

  This is a medium level function intended for use with the stellerator
  variant of neo2.
  It will set the input variable 'prop_reconstruct' to given value for
  all folders given as a list.

  WARNING: This function does not check if the given value is valid.
    This is intentionaly to avoid having to change this function if the
    range of valid values changes.

  input:
  ------
  folders: list of folders to parse.
  folder: string, determining where the folders are located. If they are
    in the local folder './' is an appropriate value.
  backup: boolean, if true the original version is saved with '~'
    appended to the name.
  value: integer, the value to use for 'prop_reconstruct'.

  output:
  -------
  No formal output.

  Side-effects:
  -------------
  Changes the files determined by input parameters, and may create
  copies of the older version.
  """
  from os.path import join

  import f90nml.parser

  parser = f90nml.parser.Parser()
  for d in folders:
    current_filename = join(folder, d.name, 'neo2.in')
    nml = parser.read(current_filename)

    print('Processing: ' + current_filename)

    if (backup):
      nml.write(current_filename + '~')

    # Adjust prop_reconstruct.
    name_value_tuples = [('propagator', 'prop_reconstruct', value)]
    change_namelist_values_for_file_object(nml, name_value_tuples)
    nml.write(current_filename, True)

def change_namelist_values_for_files(folder: str, subfolder_pattern: str, name_value_tuples: list, backup: bool):
  """Change values of namelists files in subfolders.

  \todo Allow setting of folder-dependent quantities (i.e. depending on
    flux surface).

  input:
  ------
  folder: a string, giving the location where to look for subfolders.
  subfolder_pattern: a string, giving names of subfolders where namelist
    files are located, as a regular expression, e.g 'es_*'.
  name_value_tuples: list of tuples (namelist, element, value) of values
    to set.
    Note that this is the same for all namelist files, i.e. at the
    moment it is not possible to vary values based on flux surface.
  backup: boolean, if true, a backup of the original input file is
    stored with added '~' at the end.
  """
  from os.path import join
  from pathlib import Path

  import f90nml.parser

  # Get all objects in the folder
  p = Path(folder)
  #~ # Filter for the directories
  #~ p = [x for x in p.iterdir() if x.is_dir()]
  # Filter for the name of the subdirectories
  folders = list(p.glob(subfolder_pattern))

  parser = f90nml.parser.Parser()

  for d in folders:
    current_filename = join(folder, d.name, 'neo2.in')
    nml = parser.read(current_filename)

    print('Processing: ' + current_filename)

    if (backup):
      nml.write(current_filename + '~')

    change_namelist_values_for_file_object(nml, name_value_tuples)
    nml.write(current_filename, True)

def change_namelist_values_for_file_object(fortran_namelist_file_object, name_value_tuples: list):
  """Change values of a fortran namelist object.

  Change multiple values of a fortran namelist object.

  \note This is thought of a low-level class, as a user, you probably
  want to use change_namelist_values_for_files or set_neo2in.

  input:
  ------
  fortran_namelist_file_object: a parsed fortran namelist file, as
    produced by f90nml.parser.Parser(). The values of this object are
    changed.
  name_value_tuples: list of tuples (namelist, element, value) of values
    to set.
  """
  for k in name_value_tuples:
    change_namelist_value_for_file_object(fortran_namelist_file_object, k)

def change_namelist_value_for_file_object(fortran_namelist_file_object, name_value_tuple: tuple):
  """Change one value of a fortran namelist file object.

  \note This is thought of a low-level class, as a user, you probably
  want to use change_namelist_values_for_files or set_neo2in.

  input:
  ------
  fortran_namelist_file_object: a parsed fortran namelist file, as
    produced by f90nml.parser.Parser(). The value of this object is
    changed.
  name_value_tuple: tuple (namelist, element, value) of value to set.
  """
  fortran_namelist_file_object[name_value_tuple[0]][name_value_tuple[1]] = name_value_tuple[2]

def get_list_unsucessful_runs(folder: str, subfolder_pattern: str, file_to_check: str):
  """
  \brief Return a list of unsucessful runs.

  This function will return a list of unsucessful runs, for subfolders
  of a given pattern in a given folder.
  'Unsucessful run' is thereby determined based on the absence of a file
  with given name from the subfolder. The list returned will be a list
  of the names of the subfolders.

  input:
  ------
  folder: string, with the path to the folder where to look for
    subfolders. You can use './' for the current working directory.
  subfolder_pattern: string which describes the subfolders. This may
    contain wildcards, e.g. 'es_*'.
    Attention: the pattern should not include files that are also
    present in the directory, as there is no check, to operate on
    folders only.
  file_to_check: string, name of the file which absence indicates an
    unsucessful run.

  output:
  -------
  List, containing strings with the name of the subfolders (not the full
  path) which contain runs which where not sucessfull, i.e. which do
  _not_ contain 'file_to_check'.
  If all runs were sucessful, then the list will be empty.
  """
  from os.path import join
  from pathlib import Path

  unsucessful_runs = []

  p = Path(folder)
  folders = list(p.glob(subfolder_pattern))

  for d in folders:
    current_filename = join(folder, d.name, file_to_check)

    try:
      with open(current_filename) as f:
        continue
    except FileNotFoundError:
      unsucessful_runs.append(d.name)

  unsucessful_runs.reverse() # Implementation detail - return list in increasing order.

  return unsucessful_runs


def append_list_unsucessful_runs(folder: str, subfolder_pattern: str, file_to_check: str, infilename: str, outfilename: str):
  """
  Determine list of unsucessful runs, and add list to file.

  input:
  ------
  folder: string, with the path to the folder where to look for
    subfolders. You can use './' for the current working directory.
  subfolder_pattern: string which describes the subfolders. This may
    contain wildcards, e.g. 'es_*'.
    Attention: the pattern should not include files that are also
    present in the directory, as there is no check, to operate on
    folders only.
  file_to_check: string, name of the file which absence indicates an
    unsucessful run.
  """

  unsucessful_runs = get_list_unsucessful_runs(folder, subfolder_pattern, file_to_check)

  with open(infilename,'r') as infile, open(outfilename,'w') as outfile:
    for line in infile:
      outfile.write(line)
      if len(line.split()) > 0:
        if line.split()[0] == 'Queue':
          break

    for run in unsucessful_runs:
      if run == 'es_0p00000':
        # Skip axis
        continue

      outfile.write('  {} \\\n'.format(run))


def get_runcompletion_from_output_par(outputfilename: str):
  """ Return percentage and list of propagators started.

  Helper function to determine completion of a neo-2 par (Gernots version)
  run. Gives the percentage of started propagators, and a list of the
  tags of those.

  note: uses so far started propagators. Percentage of finished
  propagators could be calculated by reducing the number of started
  propagators by the number of mpi processes-1, as each process, except
  the scheduler, has one propagator started.
  As this would require another scanning of the file, to get the number
  of mpi-processes reliable, this was considered unnecessary for now.

  input:
  ------
  outputfilename: string, the name of the file that contains the output
    which should be parsed.

  output:
  -------
  list, with first element is percantage of started propagators as
    float, and the second is a list of the propagator tags (integers).
  """

  with open(outputfilename) as f:
      lines = f.readlines()

  propagators = [match.strip('\n') for match in lines if "Level placement for propagator" in match]

  index_propagators = [int(k.split()[-1]) for k in propagators]

  nr_propagators = max(index_propagators)

  finished_propagators = [int(match.strip('\n').split()[-1]) for match in lines if "propagator tag" in match]
  # To only have unique values
  finished_propagators = list(set(finished_propagators))

  # -1 in the denumerator, because tag '1' not used(?)
  # numerator should be reduced by number of mpi-processes -1 to get
  # number of propagators finished.
  return [float(len(finished_propagators))/float(nr_propagators-1), finished_propagators]


def replace_in_file(infilename: str, outfilename, dic: dict):
  """
  Replace in file tokens with values given as dict; write to new file.

  Within a given file, replace all occurences of given tokens with
  corresponding values. Tokens and values are given as dictionary, with
  the token being the keys of the dictionary. The result is written to
  a new file.

  Example:
    dic = {
      '<S_TOKEN>': s,
      '<M_T_TOKEN>': M_t,
      '<VTH_TOKEN>': vth,
      '<EPSM_TOKEN>': epsm
      }

    replace_in_file('driftorbit.in.template', 'driftorbit.in', dict)

  input:
  ------
  infilename: string, name of the file which contains the tokens.
  outfilename: string, name of the file to which to write the result,
    i.e. the infile with the tokens replaced.
  dic: dictionary, the keys are the tokens, which are replaced with the
    corresponding values.

  output:
  -------
  none

  sideeffects:
  ------------
  Creates new file. Overwrittes file, if already a file with the given
  name exits.
  """
  import re

  pattern = re.compile('|'.join(dic.keys()))

  with open(infilename,'r') as infile, open(outfilename,'w') as outfile:
    for line in infile:
      result = pattern.sub(lambda x: str(dic[x.group()]), line)
      outfile.write(result)


if __name__ == "__main__":
  import matplotlib.pyplot as plt
  import sys

  # Check input arguments for input filename.
  if (len(sys.argv) >= 2):
    scanlogfile = sys.argv[1]
  else:
    scanlogfile = 'scan.log'

  if (len(sys.argv) >= 3):
    id_scan = sys.argv[2]
  else:
    id_scan = '1708'

  with open(scanlogfile) as f:
    lines = f.readlines()

  num_runs = 100

  memory_consumption = get_memory_consumption_of_scan(lines, id_scan, num_runs)

  plt.plot(memory_consumption)
  plt.xlabel('# Process')
  plt.ylabel('Memory in MB')
  plt.savefig('memory_consumption_over_process')
