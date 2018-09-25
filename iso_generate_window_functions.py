#!/usr/bin/env python

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
tessel_car=p['tessel_car']


auxDir = 'auxMaps_%s/'%pixel
iso_map_utils.create_directory(auxDir)


mapDir='maps_%s/'%pixel

if len(sys.argv)> 2:
    iii=int(sys.argv[2])
    mapDir='maps_%s_%03d/'%(pixel,iii)


fName='split_%d_%s.fits'%(0,freqTags[0])
template=iso_map_utils.read_map(mapDir,fName,pixel)[0]

if pixel=='healpix':
    iso_window_utils.mask_healpix(auxDir,p,template)

if pixel=='car':
    iso_window_utils.rectangle_mask_car(auxDir,p,template)

