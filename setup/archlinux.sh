#!/bin/bash

sudo -Sy --needed git cmake make ninja patch gcc gcc-fortran \
    openblas hdf5 netcdf netcdf-fortran suitesparse boost openmpi \
    fftw gsl python python-numpy debugedit fakeroot

pushd /tmp
    git clone https://aur.archlinux.org/gklib.git
    pushd gklib
        makepkg -si
    popd

    git clone https://aur.archlinux.org/metis.git
    pushd metis
        makepkg -si
    popd
popd
