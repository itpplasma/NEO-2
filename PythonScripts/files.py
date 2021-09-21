#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Michael Draxler
"""

import h5py
__version__="0.07"
import numpy as np

def compare2hdf5(path1,path2):
    ##Instance check

    b=dict()
    if isinstance(path1,(h5py.File,h5py.Group)):
        for i in path1:
            b.update(compare2hdf5(path1.get(i),path2.get(i,'Different_version')))

    elif isinstance(path1,h5py.Dataset):
        try:
            if path1.value==path2.value:
            #print('+++')

                pass
            else:
                b[path1.name]=[path1.value,path2.value]
            #print("---",[path1.name],[path1.value,path2.value])
        except ValueError:
            if np.allclose(path1,path2):
                pass
            else:
                b[path1.name]=[path1.value,path2.value]
        except AttributeError:
            if path2=='Different_version':
                pass
            else:
                raise AttributeError('file2 has no attribute value')
    else:
        raise TypeError('Input must be h5py')

    return b

def comparehdf5(*paths):
    ### Input paths which should be compared

    if len(paths)<2:
        raise TypeError('Input expected at least 2 arguments')
    diff=dict()
    init=True
    for path in paths:
        if init:
            file1=h5py.File(path)
            init=False
            continue
        file2=h5py.File(path)
        dict_new=compare2hdf5(file1,file2) # Compare always first File with the
        for key,value in dict_new.items():
            if key in diff:
                diff[key].append(value[1]) # Append only the new value
            else:
                diff[key]=[*value] ## Append both values

        file2.close()

    file1.close()
    return diff

#class h5pz(h5py):

 #   comparehdf5=comparehdf5
