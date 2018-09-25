#!/usr/bin/env python

import healpy as hp
import numpy as np
import pylab as plt
import iso_dict
import iso_ps_utils
import iso_map_utils
import iso_window_utils
import iso_mode_coupling_utils
import sys
import os



p = iso_dict.flipperDict()
p.read_from_file(sys.argv[1])

freqTags=p['freqTags']
pixel=p['pixelisation']

auxDir = 'auxMaps_%s/'%pixel
mcmDir = 'mcm_%s/'%pixel
lmax= p['lmax']
type=p['type']

iso_map_utils.create_directory(mcmDir)

binLo,binHi,binSize= iso_ps_utils.read_binning_file(p['binningFile'],p['lmax'])


winList=auxDir+'window_list.txt'
spectraList=iso_window_utils.get_spectra_list(winList)

for spec in spectraList:
    print spec
    i,f1,f2=spec[0],spec[1],spec[2]
    ell,Wl_array=iso_ps_utils.get_window_beam_array(p,f1,f2,lmax)
    iso_mode_coupling_utils.get_mode_coupling(auxDir,mcmDir,i,f1,f2,Wl_array,binLo,binHi,binSize,lmax,type,pixel)
                

