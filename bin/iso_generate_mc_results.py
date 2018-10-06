#!/usr/bin/env python

# This executable combine spectra from simulations and give you their mean and std, this allow to test if the estimator is unbiased.
# It uses the following arguments in the dictionnary file:
# 'nSplit': the number of split of data (each split will have a rms noise level of sqrt(nSplit)*rms)
# 'pixelisation': either car or Healpix, the pixellisation to use
# 'useNaMaster': you can use NaMaster instead of mosaic to compute the spectra, this allows to check the agreement between the two codes
# 'hdf5': when you do a lot of simulations, you will generate a lot of spectra, this argument solve this problem by storing the spectra in the hdf5 format
# 'type': 'Dl' or 'Cl', the type for spectra you want to reconstruct
# 'removeMean': you can remove the mean of the map before doing SPHT
# 'iStart': the starting index number of the simulations
# 'iStop': the final index number of the simulations

import healpy as hp
import numpy as np
import pylab as plt
import iso_spectra_utils
import iso_map_utils
import iso_monte_carlo_utils
import iso_dict
import sys
import h5py
import os

p = iso_dict.flipperDict()
p.read_from_file(sys.argv[1])

nSplits= p['nSplit']
pixel=p['pixelisation']
naMaster=p['useNaMaster']
hdf5=p['hdf5']
type=p['type']
removeMean=p['removeMean']

mapDir= 'maps_%s/'%pixel
auxDir = 'auxMaps_%s/'%pixel
mcmDir= 'mcm_%s/'%pixel
specDir = 'spectra_%s'%pixel
mcDir=  'montecarlo_%s'%pixel

if naMaster==True:
    specDir +='_namaster'
    mcDir +='_namaster'
if removeMean==True:
    specDir +='_mean_removed'
    mcDir +='_mean_removed'

iStart=p['iStart']
iStop=p['iStop']

winList=auxDir+'window_list.txt'

iso_monte_carlo_utils.get_mean_and_std(specDir,mcDir,iStart,iStop,winList,nSplits,type,hdf5)
