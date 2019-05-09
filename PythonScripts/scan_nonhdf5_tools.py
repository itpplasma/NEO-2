#!/usr/bin/env python3

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
