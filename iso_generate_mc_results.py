#!/usr/bin/env python

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
binningFile=p['binningFile']
type=p['type']
lmax=p['lmax']
removeMean=p['removeMean']
theoryFile=p['clfile_trunc']


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
