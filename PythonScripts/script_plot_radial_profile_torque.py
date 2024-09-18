#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 09:06:41 2019

@author: Rico Buchholz

Script to create plots of torque density and diffusion coefficients as
function of boozer s for all subfolders 'n1_vshift*', and to write the
data to a text file.
A total of five plots for each folder will be created in png and eps
format.
The plotted quantities are
- TphiNA_spec
- D11_AX
- D11_NA
- D12_AX
- D12_NA
The radial coordinate and the first four entries of the diffusion
coefficients are also written to a text file in column format.
Note that ff the run had only two species then the diffusion
coefficients have only four entries in total.
"""

if __name__ == "__main__":
  from os.path import join
  from pathlib import Path
  from neo2_util import get_hdf5file
  import matplotlib.pyplot as plt

  p = Path('./')
  folders = list(p.glob('n1_vshift*'))
  for d in folders:
    plt.clf()
    try:
      with get_hdf5file(join(d.name, 'neo2_multispecies_out.h5')) as h:

        # torque
        plt.plot(h['boozer_s'], h['TphiNA_tot'], h['boozer_s'], h['TphiNA_spec'], [0, 1], [0, 0])
        plt.title('Torque over radius')
        #~ plt.legend(['dlamdv','dlamdv0'])
        #~ plt.xlim(0,5)
        #~ plt.ylim(0,5)
        plt.xlabel('s')
        plt.ylabel('NTV total torque')
        plt.savefig('NTV_total_torque'+d.name+'.eps')
        plt.savefig('NTV_total_torque'+d.name+'.png')
        #~ plt.show()

        # D11_AX
        plt.plot(h['boozer_s'], h['D11_AX'])
        plt.title('D11_AX over radius')
        plt.xlabel('s')
        plt.ylabel('D11_AX')
        plt.savefig('D11_AX_'+d.name+'.eps')
        plt.savefig('D11_AX_'+d.name+'.png')

        # D11_NA
        plt.plot(h['boozer_s'], h['D11_NA'])
        plt.title('D11_NA over radius')
        plt.xlabel('s')
        plt.ylabel('D11_NA')
        plt.savefig('D11_NA_'+d.name+'.eps')
        plt.savefig('D11_NA_'+d.name+'.png')

        # D12_AX
        plt.plot(h['boozer_s'], h['D12_AX'])
        plt.title('D12_AX over radius')
        plt.xlabel('s')
        plt.ylabel('D12_AX')
        plt.savefig('D12_AX_'+d.name+'.eps')
        plt.savefig('D12_AX_'+d.name+'.png')

        # D12_NA
        plt.plot(h['boozer_s'], h['D12_NA'])
        plt.title('D12_NA over radius')
        plt.xlabel('s')
        plt.ylabel('D12_NA')
        plt.savefig('D12_NA_'+d.name+'.eps')
        plt.savefig('D12_NA_'+d.name+'.png')

        with open('D_'+d.name+'.dat','w') as f:
          for x in range(len(h['D12_NA'])):
            f.write('{:16.8e} {:16.8e} {:16.8e} {:16.8e} {:16.8e} {:16.8e} {:16.8e} {:16.8e} {:16.8e}\n'.format(h['boozer_s'][x], h['D11_NA'][x, 0], h['D11_NA'][x, 1], h['D11_NA'][x, 2], h['D11_NA'][x, 3], h['D12_NA'][x, 0], h['D12_NA'][x, 1], h['D12_NA'][x, 2], h['D12_NA'][x, 3]))
    except OSError:
      # Assumption: no file found for this folder - ignore it.
      pass
