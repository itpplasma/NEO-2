#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Michael Draxler
"""

import h5py
def comparehdf5(self,path1,path2):
    ##Instance check

    b=dict()
    if isinstance(path1,(h5py.File,h5py.Group)):
        for i,j in zip(path1.values(),path2.values()):
            print(i.name,b)
            b.update(comparehdf5(i,j))

    elif isinstance(path1,h5py.Dataset):
        if path1.value==path2.value:
            print('+++')
            pass
        else:
            b[path1.name]=[path1.value,path2.value]
            print("---",[path1.name],[path1.value,path2.value])
    else:
        return 1

    return b

class h5pz(h5py):

    comparehdf5=comparehdf5
