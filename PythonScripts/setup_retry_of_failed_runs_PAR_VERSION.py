#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###################################################################################################
###########################################FUNCTIONS###############################################
###################################################################################################

def get_list_failed_runs_par(folder: str, subfolder_pattern: str, file_to_check: str):
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

  failed_runs = []

  p = Path(folder)
  folders = list(p.glob(subfolder_pattern))

  for d in folders:
    current_filename = join(folder, d.name, file_to_check)

    matching_files = g.glob(current_filename)

    if matching_files:
        failed_runs.append(d.name)
    else:
        continue

  failed_runs.reverse() # Implementation detail - return list in increasing order.

  return failed_runs

###################################################################################################

def append_list_failed_runs_par(failed_runs: list, infilename: str, outfilename: str):
  """
  Determine add list of unsucessful runs to file.

  Intended for adding the list at the end of a condor_submit file.
  For this reason:
  - lines after the 'Queue' command are ignored.
  - list is written as one file per line with '\' at the end, to
    indicate line continuation.
    Note that also the last line has a '\', as this does not seem to be
    a problem for HTcondor and it makes the code simpler.

  Example:
    append_list_failed_runs(failed_runs 'submit_failed_template', 'submit_failed')

  input:
  ------
  failed_runs: list of strings, containing the names of the
    subfolders of the unsucessful runs.
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

  with open(infilename,'r') as infile, open(outfilename,'w') as outfile:
    for line in infile:
      outfile.write(line)
      if len(line.split()) > 0:
        if line.split()[0] == 'Queue':
          break

    for run in failed_runs:
      if run == 's1.000000000m03':
        # Skip axis
        continue

      outfile.write('  {} \\\n'.format(run))

###################################################################################################
###########################################MAIN####################################################
###################################################################################################

def setup_retry_of_failed_runs_PAR_VERSION(faulty_file_pattern: str = "propagator*.h5", use_record_of_failed_runs: bool = False, 
                                           eps: float = 1e-09, no_sideeffects: bool = False):
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
    use_record_of_failed_runs: bool, if True, the list of failed runs is read from record_of_failed_runs.txt
    eps: float, minimal agreement between entry in surfaces.dat and the name of the detected faulty run 
         in terms of (relative) numerical difference. The folder name is shorter compared to the flux
         label entry in sufaces.dat, which leads to a mismatch in the last digits (rounding of last digit in 
         folder name that can propagate to higher digits as well) The whole logic is first performed only with 
         strings (to avoid conversion erros). However, should the flux label of some of the faulty runs not be
         detected this way, a float number comparison of these not detected cases with the entries in surfaces.dat
         is performed. Should still no sucessfull match happen, one ought to change eps or check surfaces.dat by hand.
         (default: 1e-09)
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

    if not use_record_of_failed_runs:
        failed_runs = get_list_failed_runs_par("./", "s[1234567890]*", faulty_file_pattern)
    else:
        with open('record_of_failed_runs.txt', 'r') as file:
            failed_runs = file.read().splitlines()
            print('Using record_of_failed_runs.txt to determine failed runs.')
            for failed_run in failed_runs:
                print(failed_run)

    if len(failed_runs) == 0:
        print('No failed runs found. Nothing to do.')
        return
    
    if no_sideeffects:
        print('Only printing list of failed runs. No sideeffects.')
        for failed_run in failed_runs:
                print(failed_run)
        return
    
    # After recreating the folder in later steps, the runs are not considered as failed anymore.
    # Therefore, the list of failed runs is written to a file as backup.
    # Existing entries are not overwritten to further not loose the information.
    if not use_record_of_failed_runs:
        with open('record_of_failed_runs.txt', 'a') as file:
            for i in range(len(failed_runs)):
                file.write(failed_runs[i] + '\n')

    flux_labels_of_failed_runs = [run.split('m')[0] + '*' for run in failed_runs]
    flux_labels_of_failed_runs = [flux_label.split('s')[1] for flux_label in flux_labels_of_failed_runs]
    exponents_of_failed_runs = [int(run.split('m')[1]) for run in failed_runs]
    for i in range(len(flux_labels_of_failed_runs)):
        parts = flux_labels_of_failed_runs[i].split('.')
        flux_labels_of_failed_runs[i] = '0.' + (exponents_of_failed_runs[i]-1)*'0' + parts[0] + parts[1]

    found_surfaces_of_failed_runs = []
    with open('surfaces.dat', 'r') as surface_file:
        for line in surface_file:
            columns = line.strip().split(' ')
            flux_label_of_line = columns[0] 
            if any(fnmatch.fnmatch(flux_label_of_line, flux_label) for flux_label in flux_labels_of_failed_runs):
                found_surfaces_of_failed_runs.append(line)

    if len(found_surfaces_of_failed_runs) < len(failed_runs):
        not_matched_flux_labels = [flux_label[:-1] for flux_label in flux_labels_of_failed_runs if not any(fnmatch.fnmatch(surface.strip().split(' ')[0],flux_label) for surface in found_surfaces_of_failed_runs)]
        not_matched_flux_labels = [float(flux_label) for flux_label in not_matched_flux_labels]
        print('Up to this point found '+ str(len(found_surfaces_of_failed_runs)) + ' surfaces:')
        for surface in found_surfaces_of_failed_runs:
            print(surface, end='')
        print('Up to this point not matched ' + str(len(not_matched_flux_labels)) + ' flux_labels:')
        for flux_label in not_matched_flux_labels:
            print(flux_label)
        print('Could not find all surfaces. Convert to float to handle rounding and scientific notation format in surfaces.dat.') 
        with open('surfaces.dat', 'r') as surface_file:
            for line in surface_file:
                columns = line.strip().split(' ')
                flux_label_of_line = float(columns[0])

                if abs(flux_label_of_line) < eps:
                    reference = 1 # compare absolut difference
                else:
                    reference = flux_label_of_line # compare relative difference

                for flux_label in not_matched_flux_labels:
                    if abs(flux_label_of_line - flux_label)/reference < eps:
                        found_surfaces_of_failed_runs.append(line)
                        not_matched_flux_labels.remove(flux_label)
                        break

        if len(found_surfaces_of_failed_runs) < len(failed_runs):
            print('Could not find all surfaces, considering a relative tolerance of ' + str(eps))
            print('Possibly increase tolerance eps or check surfaces.dat for not agreeing case.')
            print('List of not matched flux_labels:')
            print(not_matched_flux_labels)
    
    if len(found_surfaces_of_failed_runs) > len(failed_runs):
        print('Data found corresponding to the not successfull runs:')
        print(found_surfaces_of_failed_runs)
        raise ValueError('Something went wrong. More surfaces found than failed runs.')

    print('Found in total ' + str(len(found_surfaces_of_failed_runs)) + ' matches between the surfaces in surfaces.dat and the ' + str(len(failed_runs)) + ' failed runs.')

    for surface in found_surfaces_of_failed_runs:
        print(surface, end='')

    with open('surfaces_failed.dat', 'w') as file:
        file.writelines(found_surfaces_of_failed_runs)

    # Make submit file for failed runs
    append_list_failed_runs_par(failed_runs,"submit_failed.template", "submit_failed")

    # Make surfaces folders for failed runs
    print('Executing create_surf.py to generate new folders of failed runs:')
    import create_surf
    create_surf.create_surfaces('surfaces_failed.dat')

    return

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Setup retry of failed runs for PAR version.')
    parser.add_argument('--faulty_file_pattern', type=str, default="propagator*.h5", help='Name of the file which PRESENCE indicates a failure.')
    parser.add_argument('--use_record_of_failed_runs', type=bool, default=False, help='If True, the list of failed runs is read from record_of_failed_runs.txt')
    parser.add_argument('--eps', type=float, default=1e-09, help='Allowed (relative) rounding difference between folder name and flux label in surfaces.dat.')
    parser.add_argument('--no_sideeffects', type=bool, default=False, help='If True, only print list of failed runs. No sideeffects.')
    args = parser.parse_args()
    setup_retry_of_failed_runs_PAR_VERSION(args.faulty_file_pattern, args.use_record_of_failed_runs, args.eps, args.no_sideeffects)
