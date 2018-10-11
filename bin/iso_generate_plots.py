#!/usr/bin/env python

# This executable plot the maps, the windows function and the spectra produced by the code.
# it uses the following arguments in the dictionnary file:
# 'nSplit': the number of split of data
# 'pixelisation': either car or Healpix, the pixellisation to use
# 'useNaMaster': you can use NaMaster instead of mosaic to compute the spectra, this allows to check the agreement between the two codes
# 'hdf5': when you do a lot of simulations, you will generate a lot of spectra, this argument solve this problem by storing the spectra in the hdf5 format
# 'type': 'Dl' or 'Cl', the type for spectra you want to reconstruct
# 'lmax': the maximum multipole we want to consider
# 'theoryFile': the theoretical lensed power spectrum from CAMB you want to compare the data to
# 'survey_mask_coordinate': for car you should specify the patch coordinates you want to analyse
# 'tessel_healpix': for healpix, you can choose to do a multiple patches analysis tesseling the healpix map with patch center the pixel center corresponding to the 'patch_nside' argument and radius 'radius', the 'cut_threshold' argument allow to get rid of patches falling in the galactic mask.
# 'freqTags': the different frequency band you wish to analyze.
# 'removeMean': you can remove the mean of the map before doing SPHT
# 'compareMosaicNaMaster': if True and if you have ran the code once with mosaic and once with namaster it will overplot the two spectra
# 'colorRange': the colorRange used to plot the map it's an array of the form [T,Q,U]
# 'useMcErrors': wether you use MC errorsbar or not


import healpy as hp
import numpy as np
import pylab as plt
import iso_dict
import iso_map_plot_utils
import iso_spectra_plot_utils
import iso_map_utils
import sys
import h5py
import os

p = iso_dict.flipperDict()
p.read_from_file(sys.argv[1])

freqTags=p['freqTags']
nSplits= p['nSplit']
pixel=p['pixelisation']
hdf5=p['hdf5']
type=p['type']
lmax=p['lmax']
theoryFile=p['clfile_trunc']
survey_mask_coordinates=p['survey_mask_coordinate']
tessel_healpix=p['tessel_healpix']
removeMean=p['removeMean']
compareMosaicNaMaster=p['compareMosaicNaMaster']
colorRange=p['colorRange']
naMaster=p['useNaMaster']
mask=p['mask']
noise=p['noise']

mapDir= 'maps_%s/'%pixel
auxDir = 'auxMaps_%s/'%pixel
mcmDir= 'mcm_%s/'%pixel
specDir = 'spectra_%s'%pixel
plotDir= 'plots_%s'%pixel
mcDir=  'montecarlo_%s'%pixel

if naMaster==True:
    specDir +='_namaster'
    mcDir+='_namaster'
    plotDir+='_namaster'
if removeMean==True:
    specDir +='_mean_removed'
    mcDir +='_mean_removed'
    plotDir +='_mean_removed'
if len(sys.argv)> 2:
    mapDir= 'maps_%s_%s/'%(pixel,sys.argv[2])
    specDir += '_%s'%(sys.argv[2])
    mcDir += '_%s'%(sys.argv[2])
    plotDir += '_%s'%(sys.argv[2])

iso_map_utils.create_directory(plotDir)

winList=auxDir+'window_list.txt'

white_noise_level={}
for f1 in freqTags:
    for f2 in freqTags:
        if f1==f2:
            white_noise_level[f1,f2,'noiseInfo']=nSplits,noise['rms_%s_T'%f1],noise['rms_%s_pol'%f2]
        else :
            white_noise_level[f1,f2,'noiseInfo']=nSplits,0,0
        
        white_noise_level[f1,f2,'beamName']=p['beam_%s_T'%f1],p['beam_%s_pol'%f2]

if p['useMcErrors'] ==False:
    mcDir=None


iso_map_plot_utils.plot_survey_map(auxDir,mapDir,plotDir,pixel,winList,mask,freqTags,survey_mask_coordinates=survey_mask_coordinates,tessel_healpix=tessel_healpix,color_range=None)
iso_map_plot_utils.plot_all_windows(auxDir,plotDir,pixel,winList,tessel_healpix=tessel_healpix)
iso_map_plot_utils.plot_all_maps(auxDir,mapDir,plotDir,pixel,winList,nSplits,color_range=colorRange,survey_mask_coordinates=survey_mask_coordinates,tessel_healpix=tessel_healpix)
iso_spectra_plot_utils.plot_all_spectra(mcmDir,specDir,plotDir,winList,nSplits,hdf5,type,theoryFile,lmax,compare_mosaic_namaster=compareMosaicNaMaster,white_noise_level=white_noise_level,mcDir=mcDir)


