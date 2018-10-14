#!/usr/bin/env python

import healpy as hp
import numpy as np
import iso_map_utils
import iso_ps_utils
import os
import iso_dict
import sys
from mpi4py import MPI
import pylab as plt

p = iso_dict.flipperDict()
p.read_from_file(sys.argv[1])

freqTags=p['freqTags']
pixel=p['pixelisation']
pixWin=p['pixWin']
dataDir=p['dataDir']

mapDir='maps_%s/'%(pixel)
iso_map_utils.create_directory(mapDir)

for freq in freqTags:
    mapFiles=p['dataMap_%s'%freq]
    for count,m in enumerate(mapFiles):
        m=iso_map_utils.read_map(dataDir,m,pixel)
        
        if pixel=='car' and pixWin==True:
            from enlib import enmap
            m= enmap.apply_window(m, pow=-1.0)

        fName='split_%d_%s.fits'%(count,freq)
        iso_map_utils.write_map(m,mapDir,fName,pixel)

