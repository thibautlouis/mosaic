#!/usr/bin/env python

# This executable compute all the power spectra of a given set of maps.
# The executable uses the following arguments in the dictionnary file:
# 'nSplit': the number of split of data
# 'pixelisation': either car or Healpix, the pixellisation to use
# 'useNaMaster': you can use NaMaster instead of mosaic to compute the spectra, this allows to check the agreement between the two codes
# 'hdf5': when you do a lot of simulations, you will generate a lot of spectra, this argument solve this problem by storing the spectra in the hdf5 format
# 'binningFile': a binningFile with three colums: binMin, binMax, binMean
# 'type': 'Dl' or 'Cl', the type for spectra you want to reconstruct
# 'lmax': the maximum multipole we want to consider
# 'survey_mask_coordinate': for car you should specify the patch coordinates you want to analyse
# 'removeMean': you can remove the mean of the map before doing SPHT
# 'thetaCut': for healpix you can use a cut in the theta integral for doing SPHT, this is particularly relevant if you want to analyze small patches
# 'niter': number of iteration for doing the SPHT, this increase their accuracy, default is 3, but if you want many spectra to generate a covariance matrix you probably good to go with 0 iterations.


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
thetaCut=p['thetaCut']
niter=p['niter']
pixWin=p['pixWin']

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
    iso_spectra_utils.get_spectra(mapDir,auxDir,mcmDir,specDir,winList,nSplits,niter,lmax,binningFile,type,hdf5,pixel,pixWin,survey_mask_coordinates=survey_mask_coordinates,removeMean=removeMean,thetaCut=thetaCut)
else:
    iso_spectra_utils.get_spectra_namaster(p,mapDir,auxDir,mcmDir,specDir,winList,nSplits,lmax,type,pixel,pixWin,nlb=50)

