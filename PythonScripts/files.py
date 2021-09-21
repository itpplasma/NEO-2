#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Michael Draxler
"""

import h5py
__version__="0.03"

def comparehdf5(path1,path2):
    ##Instance check

    b=dict()
    if isinstance(path1,(h5py.File,h5py.Group)):
        for i,j in zip(path1.values(),path2.values()):
            #print(i.name,b)
            b.update(comparehdf5(i,j))

    elif isinstance(path1,h5py.Dataset):
        if path1.value==path2.value:
            #print('+++')
            pass
        else:
            b[path1.name]=[path1.value,path2.value]
            #print("---",[path1.name],[path1.value,path2.value])
    else:
        raise TypeError('Input must be h5py')

    return b


def compare(*files):

    if len(files)<2:
        raise TypeError('Input expected at least 2 arguments')
    diff=dict()
    init=True
    for path in files:
        if init:
            path1=h5py.File(path)
            init=False
            continue
        path2=h5py.File(path)
        dict_new=comparehdf5(path1,path2)
        for key,value in dict_new.items():
            if key in diff:
                diff[key].extend(value)
            else:
                diff[key]=[*value]

        path2.close()

    path1.close()
    return diff

#class h5pz(h5py):

 #   comparehdf5=comparehdf5
