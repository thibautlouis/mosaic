#!/usr/bin/env python

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

nSplits= p['nSplit']
pixel=p['pixelisation']
naMaster=p['useNaMaster']
hdf5=p['hdf5']
type=p['type']
lmax=p['lmax']
theoryFile=p['clfile_trunc']
survey_mask_coordinates=p['survey_mask_coordinate']
tessel_healpix=p['tessel_healpix']
freqTags=p['freqTags']
removeMean=p['removeMean']
compareMosaicNaMaster=p['compareMosaicNaMaster']

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
            white_noise_level[f1,f2,'noiseInfo']=nSplits,p['rmsT_%s'%f1],p['rmsP_%s'%f2]
        else :
            white_noise_level[f1,f2]=0,0
        
        white_noise_level[f1,f2,'beamName']=p['beam_%s_T'%f1],p['beam_%s_Pol'%f2]

mcDir=None
iso_map_plot_utils.plot_survey_map(p,auxDir,mapDir,plotDir,pixel,winList,survey_mask_coordinates=survey_mask_coordinates,tessel_healpix=tessel_healpix,color_range=None)
iso_map_plot_utils.plot_all_windows(auxDir,plotDir,pixel,winList,tessel_healpix=tessel_healpix)
color_range=[400,50,50]
iso_map_plot_utils.plot_all_maps(auxDir,mapDir,plotDir,pixel,winList,nSplits,color_range=color_range,survey_mask_coordinates=survey_mask_coordinates,tessel_healpix=tessel_healpix)
#iso_spectra_plot_utils.plot_all_spectra(mcmDir,specDir,hdf5,winList,nSplits,tessel_healpix=tessel_healpix)

iso_spectra_plot_utils.plot_all_spectra(mcmDir,specDir,plotDir,winList,nSplits,hdf5,type,theoryFile,lmax,compare_mosaic_namaster=compareMosaicNaMaster,white_noise_level=white_noise_level,mcDir=mcDir)
