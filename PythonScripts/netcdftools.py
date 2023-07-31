#!/usr/bin/env python3

def plot_quantity(file_name:str, field_name:str, mode_number: int = 0):
  """Plot from a given file a given quantity.

  Plot for a given file, a specific field over its first dimension, with
  given index for second dimension (zero per default).
  For x-axis the index of the first dimension will be used.

  input:
  ------
  file_name: string, path+name of the file from which to plot.
  field_name: string, name of the field which to plot.
  mode_number: integer, optional [default=0], the second dimension index
    of the quantitiy to plot.

  output:
  -------
  No direct output.

  sideeffects:
  ------------
  Not sure if the plotting counts, as it opens a window that blocks
  further execution.
  """

  import matplotlib.pyplot as plt
  import numpy as np
  import scipy.io.netcdf as ncdf

  with ncdf.netcdf_file(file_name) as netcdf_file:
    netcdf_data = netcdf_file.variables
    data =  np.copy(netcdf_data[field_name].data)

  plt.figure(1)
  plt.clf()
  plt.title(field_name)
  plt.plot(data[:,mode_number])
  plt.show()


if __name__ == "__main__":
  pass
