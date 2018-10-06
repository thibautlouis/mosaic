#!/usr/bin/env python

# This executable generates window functions that will define the part of the map in which the power spectrum is computed
# It allows for doing the analysis on the full observed map or a set of different patches
# The window function are written in auxDir, the script also write the list of window function we should consider in auxDir.
# The executable uses the following arguments in the dictionnary file:
# 'pixelisation': either car or Healpix, the pixellisation to use
# 'freqTags': the different frequency band you wish to simulate.
# 'mask': a dictionnary, you can choose to use galactic mask and pts sources mask in T and pol, if you choose to use these mask, you need to specify their location in the mask dictionnary object defined in global.dict
# 'apo_type': the type of apodisation you want to use, you can use the standard C1,C2 apodisation.
# C1_namaster,C2_namaster and smooth_namaster are also available for healpix, in that case, the code will use NaMaster to generate the apod
# C3 is available for CAR, it's a function that was used for flat sky actpol analysis.
# 'apo_radius': this is dictionnary taking as argument the apodisation radius in degree, you can choose different values for the survey (this include the real survey window and the galactic mask, and for the point sources)
# 'tessel_healpix': for healpix, you can choose to do a multiple patches analysis tesseling the healpix map with patch center the pixel center corresponding to the 'patch_nside' argument and radius 'radius', the 'cut_threshold' argument allow to get rid of patches falling in the galactic mask.
# 'survey_mask_coordinate': for car you should specify the patch coordinates you want to analyse


import healpy as hp
import numpy as np
import pylab as plt
import iso_dict
import iso_map_utils
import iso_window_utils
import iso_apodization_utils
import sys
import os

p = iso_dict.flipperDict()
p.read_from_file(sys.argv[1])

freqTags=p['freqTags']
pixel=p['pixelisation']
mask=p['mask']
tessel_healpix=p['tessel_healpix']
survey_mask_coordinate=p['survey_mask_coordinate']
apo_type=p['apo_type']
apo_radius=p['apo_radius']


auxDir = 'auxMaps_%s/'%pixel
iso_map_utils.create_directory(auxDir)

mapDir='maps_%s/'%pixel

if len(sys.argv)> 2:
    iii=int(sys.argv[2])
    mapDir='maps_%s_%03d/'%(pixel,iii)

fName='split_%d_%s.fits'%(0,freqTags[0])
template=iso_map_utils.read_map(mapDir,fName,pixel)[0]

if pixel=='healpix':
    iso_window_utils.mask_healpix(auxDir, freqTags, pixel, tessel_healpix, mask, apo_type, apo_radius, template, plot=False)

if pixel=='car':
    iso_window_utils.rectangle_mask_car(auxDir, freqTags, pixel, survey_mask_coordinate, mask, apo_type, apo_radius, template, plot=False)

