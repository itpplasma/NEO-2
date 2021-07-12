#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 01 13:30:41 2021

@author: Rico Buchholz

Script/functions to detmine if two files are equal.

Functions can be used, or this file can be invoked as a script. If it
is invoked as a script, it will expect four commandline arguments: two
filenames, absolute accuracy and relative accuracy. If less than four
arguments are provided, a help message is provided, but no comparisons
are done.
"""

def line_to_float(line: str):
  """
  Try to convert a line of text to floating point numbers.

  This function tries to convert the given line of text to a list of
  floating point numbers. This is done by first splitting the line and
  then try to convert each element of the resulting list. Excepetions
  are not handled by this function, i.e. for unknown input you should at
  least expect ValueError to be possible.

  input:
  ------
  line: string with is considered to be a line of text, which (may)
    contain floating point numbers.

  output:
  -------
  A list of floating point numbers. Might be empty if there are no
  numbers in the string. If one of the conversions results in an
  exception then you will get the exception, it is not handled here.
  """
  ls = line.split()
  floats = []
  for l in ls:
    floats.append(float(l))

  return floats

def compare_numbers(numbers_1: list, numbers_2: list, abs_accuracy: float, rel_accuracy: float):
  """
  Compare two lists of lists of numbers.

  Check if two lists of lists of numbers are equal.
  Two lists of numbers are considered equal if they have the same number
  of elements, and if each of the numbers has the same value.
  Two numbers are considered to have the same value here if:
  1. theyre absolute difference is smaller than a given threshold.
  2. theyre difference relative to the largest value in the list is
  smaller than a given threshold.

  input:
  ------
  numbers_1, numbers_2: lists, to be compared.
  abs_accuracy: absolute accuracy to use, i.e. two numbers are
    considered to differ if |1 - n2| > abs_accuracy.
  rel_accuracy: relative accurracy to use, i.e. two numbers are
    considered to differ if |(n1 - n2)/max(n1| > rel_accuracy.

  output:
  -------
  True if the lists are the same, false if not.
  """

  if len(numbers_1) == len(numbers_2):
    for l1,l2 in zip(numbers_1, numbers_2):
      if len(l1) > 0:
        max_ = max(l1)
        if max_ == 0:
          max_ = 1.0
      else:
        continue

      if len(l1) == len(l2):
        for f1,f2 in zip(l1, l2):
          if (abs(f1 - f2) > abs_accuracy or abs((f1 - f2)/max_) > rel_accuracy):
            return False
      else:
        return False
  else:
    return False

  return True

def compare_number_files(file_1: str, file_2: str, abs_accuracy: float, rel_accuracy: float):
  """
  Check if two files containing numbers are the same.

  Check if two files containing numbers are the same. This means that
  all the lines with only numbers are the same. If a line contains data
  that can not be converted into a number, then the line is ignored.
  For equality same principles as for compare_numbers apply.

  input:
  ------
  file_1, file_2: strings containing file name and path of the files to
    compare.
  abs_accuracy: absolute accuracy to use, i.e. two numbers are
    considered to differ if |1 - n2| > abs_accuracy.
  rel_accuracy: relative accurracy to use, i.e. two numbers are
    considered to differ if |(n1 - n2)/max(n1| > rel_accuracy.

  output:
  -------
  True if the files have the same data, false otherwise.
  """
  # Read in the files.
  try:
    with open(file_1) as f:
      lines_1 = f.readlines()
  except FileNotFoundError:
    print('File: "' + file_1 + '" could not be loaded')
    return False

  try:
    with open(file_2) as f:
      lines_2 = f.readlines()
  except FileNotFoundError:
    print('File: "' + file_2 + '" could not be loaded')
    return False

  numbers_1 = []
  for l in lines_1:
    try:
      numbers_1.append(line_to_float(l))
    except ValueError:
      pass

  numbers_2 = []
  for l in lines_2:
    try:
      numbers_2.append(line_to_float(l))
    except ValueError:
      pass

  return compare_numbers(numbers_1, numbers_2, abs_accuracy, rel_accuracy)

def compare_text_number_files(file_1: str, file_2: str, abs_accuracy: float, rel_accuracy: float):
  """
  Check if two files containing numbers and text are the same.

  Check if two files, which may contain numbers and text, are the same.
  This is done on a line-by-line basis: try to convert line into numbers
  if this is possible, compare them, if not compare the text.

  input:
  ------
  file_1, file_2: strings, names and paths of the files which to
    compare.
  abs_accuracy: absolute accuracy to use, i.e. two numbers are
    considered to differ if |1 - n2| > abs_accuracy.
  rel_accuracy: relative accurracy to use, i.e. two numbers are
    considered to differ if |(n1 - n2)/max(n1| > rel_accuracy.

  output:
  -------
  True if the files are equal, false if not.
  """
  # Read in the files.
  try:
    with open(file_1) as f:
      lines_1 = f.readlines()
  except FileNotFoundError:
    print('File: "' + file_1 + '" could not be loaded')
    return False

  try:
    with open(file_2) as f:
      lines_2 = f.readlines()
  except FileNotFoundError:
    print('File: "' + file_2 + '" could not be loaded')
    return False

  numbers_1 = []
  text_1 = []
  for l in lines_1:
    try:
      numbers_1.append(line_to_float(l))
    except ValueError:
      text_1.append(l)

  numbers_2 = []
  for l in lines_2:
    try:
      numbers_2.append(line_to_float(l))
    except ValueError:
      text_2.append(l)

  if len(text_1) == len(text_2):
    for l1, l2 in zip(text_1, text_2):
      if l1 != l2:
        return False
  else:
    return False

  return compare_numbers(numbers_1, numbers_2, abs_accuracy, rel_accuracy)

def compare_text_files(file_1: str, file_2: str, abs_accuracy: float, rel_accuracy: float):
  """
  check if two text files are equal.

  Check if two text files are equal. They are considered equal if they
  have the same amount of lines and each line is the same.

  input:
  ------
  file_1, file_2: strings, names and paths of the files which to
    compare.
  abs_accurac, rel_accuracy: floats, not used. Present to have same
  interface as other comparison routines.

  output:
  -------
  True if the files are equal, false otherwise.
  """

  if len(lines_1) == len(lines_2):
    for k,j in zip(lines_1, lines_2):
      if k != j:
        return False

    return True
  else:
    return False

def compare_data_files(file_1: str, file_2: str, abs_accuracy: float, rel_accuracy: float, switch_comparison_type: int):
  """
  Compare content of two data files.

  input:
  ------
  file_1: string, name of first file.
  file_2: string, name of second file.
  abs_accuracy: floating point number, absolute accuracy to use for equality.
  rel_accuracy: floating point number, relative accuracy to use for equality.
  switch_comparison_type: integer switch.
    0 - do not compare strings
    1 - do compare strings
    2 - compare as strings

  output:
  -------
  True if files are considered equal, false if not.
  """

  tel = {0: compare_number_files, 1: compare_text_number_files, 2: compare_text_files}

  return tel[switch_comparison_type](file_1, file_2, abs_accuracy, rel_accuracy)

if __name__ == "__main__":
  import sys

  if (len(sys.argv) < 4):
    print("Usage:")
    print("./compare_data_files.py file_1 file_2 abs_accuracy rel_accuracy")

    sys.exit(1)
  else:
    file_1 = sys.argv[1]
    file_2 = sys.argv[2]
    abs_accuracy = float(sys.argv[3])
    rel_accuracy = float(sys.argv[4])

    if compare_data_files(file_1, file_2, abs_accuracy, rel_accuracy, 0):
      sys.exit(0)
    else:
      sys.exit(1)
