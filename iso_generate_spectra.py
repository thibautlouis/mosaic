#!/usr/bin/env python

import healpy as hp
import numpy as np
import pylab as plt
import iso_spectra_utils
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
survey_mask_coordinates=p['survey_mask_coordinate']
removeMean=p['removeMean']
theta_cut=p['theta_cut']
niter=p['niter']

mapDir= 'maps_%s/'%pixel
auxDir = 'auxMaps_%s/'%pixel
mcmDir= 'mcm_%s/'%pixel
specDir = 'spectra_%s'%pixel

if naMaster==True:
    specDir +='_namaster'
if removeMean==True:
    specDir +='_mean_removed'
if len(sys.argv)> 2:
    mapDir= 'maps_%s_%s/'%(pixel,sys.argv[2])
    specDir += '_%s'%(sys.argv[2])


winList=auxDir+'window_list.txt'

if naMaster==False:
    iso_spectra_utils.get_spectra(mapDir,auxDir,mcmDir,specDir,winList,nSplits,niter,lmax,binningFile,type,hdf5,pixel,survey_mask_coordinates=survey_mask_coordinates,removeMean=removeMean,theta_cut=theta_cut)
else:
    iso_spectra_utils.get_spectra_namaster(p,mapDir,auxDir,mcmDir,specDir,winList,nSplits,lmax,type,pixel,nlb=50)

