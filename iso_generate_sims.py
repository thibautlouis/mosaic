#!/usr/bin/env python

import healpy as hp
import numpy as np
import iso_map_utils
import iso_ps_utils
import iso_map_plot_utils
import os
import iso_dict
import sys
from mpi4py import MPI
import pylab as plt


p = iso_dict.flipperDict()
p.read_from_file(sys.argv[1])

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

iStart = p['iStart']
iStop = p['iStop']

delta = (iStop - iStart)/size

if delta == 0:
    raise ValueError, 'Too many processors for too small a  loop!'

iMin = iStart+rank*delta
iMax = iStart+(rank+1)*delta

if iMax>iStop:
    iMax = iStop
elif (iMax > (iStop - delta)) and iMax <iStop:
    iMax = iStop

clfile=p['clfile'] # lensed cl from CAMB
freqTags=p['freqTags']
nSplit= p['nSplit']
pixel=p['pixelisation']

if pixel=='healpix':
    nside= p['nside']
    pixArea= hp.pixelfunc.nside2pixarea(nside,degrees=True)
    lth, cl_TT_th, cl_EE_th, cl_BB_th, cl_TE_th = iso_ps_utils.read_clth_file(clfile)
if pixel=='car':
    from enlib import enmap,powspec,curvedsky,enplot
    ncomp=3
    ps=powspec.read_spectrum(clfile)[:ncomp,:ncomp]

for iii in xrange(iMin,iMax):
    
    mapDir='maps_%s_%03d/'%(pixel,iii)
    iso_map_utils.create_directory(mapDir)
    
    if pixel=='healpix':
        alms=hp.sphtfunc.synalm((cl_TT_th, cl_EE_th, cl_BB_th, cl_TE_th),new=True)
    if pixel=='car':
        alms=curvedsky.rand_alm(ps)

    for f in freqTags:
        alm_beamed=alms.copy()
        
        l,fl_T=np.loadtxt(p['beam_%s_T'%f],unpack=True)
        l,fl_Pol=np.loadtxt(p['beam_%s_Pol'%f],unpack=True)

        alm_beamed[0]=hp.sphtfunc.almxfl(alms[0], fl_T)
        alm_beamed[1]=hp.sphtfunc.almxfl(alms[1], fl_Pol)
        alm_beamed[2]=hp.sphtfunc.almxfl(alms[2], fl_Pol)
        
        
        if pixel=='healpix':
            m=hp.sphtfunc.alm2map(alm_beamed, nside)
        if pixel=='car':
            template= enmap.read_map('%s'%p['template'])
            m=curvedsky.alm2map(alm_beamed, template)

        for s in range(nSplit):
            maps=iso_map_utils.add_noise(m,p['rmsT_%s'%f],p['rmsP_%s'%f],nSplit,pixel)
            
            fName='split_%d_%s.fits'%(s,f)
            
            iso_map_utils.write_map(maps,mapDir,fName,pixel)


