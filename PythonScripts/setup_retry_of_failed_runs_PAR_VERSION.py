#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###################################################################################################
###########################################FUNCTIONS###############################################
###################################################################################################

def get_list_unsucessful_runs_par(folder: str, subfolder_pattern: str, file_to_check: str):
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
  file_to_check: string, name of the file which PRESENCE indicates an
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
  import glob as g

  unsucessful_runs = []

  p = Path(folder)
  folders = list(p.glob(subfolder_pattern))

  for d in folders:
    current_filename = join(folder, d.name, file_to_check)

    matching_files = g.glob(current_filename)

    if matching_files:
        unsucessful_runs.append(d.name)
    else:
        continue

  unsucessful_runs.reverse() # Implementation detail - return list in increasing order.

  return unsucessful_runs

###################################################################################################

def append_list_unsucessful_runs_par(folder: str, subfolder_pattern: str, file_to_check: str, infilename: str, outfilename: str):
  """
  Determine list of unsucessful runs, and add list to file.

  Intended for adding the list at the end of a condor_submit file.
  For this reason:
  - lines after the 'Queue' command are ignored.
  - list is written as one file per line with '\' at the end, to
    indicate line continuation.
    Note that also the last line has a '\', as this does not seem to be
    a problem for HTcondor and it makes the code simpler.

  Example:
    append_list_unsucessful_runs('./', 'es_*', 'neo2_multispecies_out.h5', 'submit_failed_template', 'submit_failed')

  input:
  ------
  folder: string, with the path to the folder where to look for
    subfolders. You can use './' for the current working directory.
  subfolder_pattern: string which describes the subfolders. This may
    contain wildcards, e.g. 'es_*'.
    Attention: the pattern should not include files that are also
    present in the directory, as there is no check, to operate on
    folders only.
  file_to_check: string, name of the file which PRESENCE indicates an
    unsucessful run.
  infilename: string, name (and path) of file to which to append the
    list. Lines after a line with 'Queue' a the begining (whitespace
    allowed?), are ignored.
  outfilename: string, name (and path) of file to which to write the
    result.

  output:
  -------
  none

  sideeffects:
  ------------
  Creates new file, overwritten if it already exists.
  """

  unsucessful_runs = get_list_unsucessful_runs_par(folder, subfolder_pattern, file_to_check)

  with open(infilename,'r') as infile, open(outfilename,'w') as outfile:
    for line in infile:
      outfile.write(line)
      if len(line.split()) > 0:
        if line.split()[0] == 'Queue':
          break

    for run in unsucessful_runs:
      if run == 's1.000000000m03':
        # Skip axis
        continue

      outfile.write('  {} \\\n'.format(run))

###################################################################################################
###########################################MAIN####################################################
###################################################################################################

def setup_retry_of_failed_runs_PAR_VERSION(faulty_file_pattern: str = "propagator*.h5", use_failed_record: bool = False, 
                                           eps: float = 1e-10, no_sideeffects: bool = False):
    """

    Make new submit file and reset folders for failed surfaces to retry the runs.

    Presence of propagator*.h5 is usually the faulty pattern as in a successful 
    run cylce, these files are deleted during clean up (reconstruct = 3)
    In contrast to the QL version, the PAR version does not have a specific
    output file missing in case of a failed run, but they are incomplete/deprecated.
    It is therefore easier to detect the presence of files that should not exist in case of success.

    example: setup_retry_of_failed_runs_PAR_VERSION(faulty_file_pattern = "propagator*.h5", minimum_label_length = 11)
    or in terminal with default values: python3 setup_retry_of_failed_runs_PAR_VERSION.py

    input:
    ------
    faulty_file_pattern: string, name of the file which PRESENCE indicates a failure.
    use_failed_record: bool, if True, the list of failed runs is read from record_of_failed_runs.txt
    eps: float, minimal agreement between entry in surfaces.dat and the name of the detected faulty run 
         in terms of (relative) numerical difference. The folder name is shorter compared to the flux
         label entry in sufaces.dat, which leads to a mismatch in the last digits (rounding of last digit in 
         folder name that can propagate to higher digits as well) The whole logic is first performed only with 
         strings (to avoid conversion erros). However, should the flux label of some of the faulty runs not be
         detected this way, a float number comparison of these not detected cases with the entries in surfaces.dat
         is performed. Should still no sucessfull match happen, one ought to change eps or check surface.dat by hand.
         (default: 1e-10)
    no_sideeffects: bool, if True, only print list of failed runs. No sideeffects.
    
    output:
    -------
    none

    sideeffects:
    ------------
    Creates new submit file submit_failed from submit_failed.template
    Creates new surface file surfaces_failed.dat from data extracted from surfaces.dat
    Uses surfaces_failed.dat to create new surfaces folders for failed runs with create_surf.py
    """

    # Make surfaces.dat file for failed runs
    import fnmatch

    if not use_failed_record:
        list_unsucessful_runs = get_list_unsucessful_runs_par("./", "s[1234567890]*", faulty_file_pattern)
    else:
        with open('record_of_failed_runs.txt', 'r') as file:
            list_unsucessful_runs = file.read().splitlines()
            print('Using record_of_failed_runs.txt to determine failed runs.')
            print(list_unsucessful_runs)

    if len(list_unsucessful_runs) == 0:
        print('No failed runs found. Nothing to do.')
        return
    
    if no_sideeffects:
        print('Only printing list of failed runs. No sideeffects.')
        print(list_unsucessful_runs)
        return
    
    # After recreating the folder in later steps, the runs are not considered as failed anymore.
    # Therefore, the list of failed runs is written to a file as backup.
    # Existing entries are not overwritten to further not loose the information.
    if not use_failed_record:
        with open('record_of_failed_runs.txt', 'a') as file:
            for i in range(len(list_unsucessful_runs)):
                file.write(list_unsucessful_runs[i] + '\n')

    list_unsucessful_flux_labels = [run.split('m')[0] + '*' for run in list_unsucessful_runs]
    list_unsucessful_flux_labels = [run.split('s')[1] for run in list_unsucessful_flux_labels]
    for i in range(len(list_unsucessful_flux_labels)):
        parts = list_unsucessful_flux_labels[i].split('.')
        list_unsucessful_flux_labels[i] = '0.'+parts[0]+parts[1]

    data_to_extract = []
    with open('surfaces.dat', 'r') as surface_file:
        for line in surface_file:
            columns = line.strip().split(' ')
            flux_label = columns[0] 
            if any(fnmatch.fnmatch(flux_label, run) for run in list_unsucessful_flux_labels):
                data_to_extract.append(line)

    if len(data_to_extract) < len(list_unsucessful_flux_labels):
        list_not_found = [run[:-1] for run in list_unsucessful_flux_labels if run not in [line.split(' ')[0] for line in data_to_extract]]
        print('Could not find all surfaces. Convert to float to handle rounding.')
        list_not_found = [float(run) for run in list_not_found]
        with open('surfaces.dat', 'r') as surface_file:
            for line in surface_file:
                columns = line.strip().split(' ')
                flux_label = float(columns[0])

                if abs(flux_label) < eps:
                    reference = 1 # compare absolut difference
                else:
                    reference = flux_label # compare relative difference

                if any(abs(flux_label - run)/reference < eps for run in list_not_found):
                    data_to_extract.append(line)

        if len(data_to_extract) < len(list_unsucessful_flux_labels):
            print('Could not find all surfaces, considering a relative tolerance of .' + eps)
            print('Possibly increase tolerance eps or check surface.dat for not agreeing case.')
            print('List of not found surfaces:')
            print(list_not_found)
    
    if len(data_to_extract) > len(list_unsucessful_flux_labels):
        raise ValueError('Something went wrong. More surfaces found than failed runs.')

    print('Found ' + str(len(data_to_extract)) + ' surfaces in surfaces.dat for failed runs.')

    for item in data_to_extract:
        print(item, end='')

    with open('surfaces_failed.dat', 'w') as file:
        file.writelines(data_to_extract)

    # Make submit file for failed runs
    append_list_unsucessful_runs_par("./", "s[1234567890]*", faulty_file_pattern,"submit_failed.template", "submit_failed")

    # Make surfaces folders for failed runs
    import create_surf
    create_surf.create_surfaces('surfaces_failed.dat')

    return

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Setup retry of failed runs for PAR version.')
    parser.add_argument('--faulty_file_pattern', type=str, default="propagator*.h5", help='Name of the file which PRESENCE indicates a failure.')
    parser.add_argument('--use_failed_record', type=bool, default=False, help='If True, the list of failed runs is read from record_of_failed_runs.txt')
    parser.add_argument('--eps', type=float, default=1e-10, help='Allowed (relative) rounding difference between folder name and flux label in surfaces.dat.')
    parser.add_argument('--no_sideeffects', type=bool, default=False, help='If True, only print list of failed runs. No sideeffects.')
    args = parser.parse_args()
    setup_retry_of_failed_runs_PAR_VERSION(args.faulty_file_pattern, args.use_failed_record, args.eps, args.no_sideeffects)