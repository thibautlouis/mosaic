#!/usr/bin/env python

# This executable generates CMB simulations on CAR or Healpix pixellisation from a theory power spectrum.
# It has a MPI loop in on the number of sims and uses open-mp for doing spherical harmonics transform.
# The executable uses the following arguments in the dictionnary file:
# 'clfile': the theoretical lensed power spectrum from CAMB
# 'iStart': the starting index number of the simulations
# 'iStop': the final index number of the simulations
# (simulations will be generated from iStart to iStop)
# 'pixelisation': either car or Healpix, the pixellisation to use
# 'freqTags': the different frequency band you wish to simulate.
# 'beam_'freq'_T: the beam to be applied to the temperature simulation
# 'beam_'freq'_pol: the beam to be applied to the polarisation simulation
# 'rms_'freq'_T: the white noise rms in temperature
# 'rms_'freq'_pol: the white noise rms in temperature
# 'nSplit': the number of split of data (each split will have a rms noise level of sqrt(nSplit)*rms)
#  The code create iStop-iStart folders and nSplit maps in each of them.

import healpy as hp
import numpy as np
import iso_map_utils
import iso_ps_utils
import os
import iso_dict
import sys
from mpi4py import MPI
import pylab as plt
from pixell import enmap,powspec,curvedsky,enplot

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

clfile=p['clfile']
freqTags=p['freqTags']
nSplit= p['nSplit']
pixel=p['pixelisation']
noise=p['noise']
pixWin=p['pixWin']


if pixel=='healpix':
    nside= p['nside']
    pixArea= hp.pixelfunc.nside2pixarea(nside,degrees=True)
    lth, cl_th = iso_ps_utils.read_clth_file(clfile,type='Cl')
if pixel=='car':
    ncomp=3
    ps=powspec.read_spectrum(clfile)[:ncomp,:ncomp]

for iii in xrange(iMin,iMax):
    
    mapDir='maps_%s_%03d/'%(pixel,iii)
    iso_map_utils.create_directory(mapDir)
    
    if pixel=='healpix':
        alms=hp.sphtfunc.synalm((cl_th['TT'], cl_th['EE'], cl_th['BB'], cl_th['TE']),new=True)
    if pixel=='car':
        alms=curvedsky.rand_alm(ps)

    for f in freqTags:
        
        alm_beamed=alms.copy()
        
        l,fl_T=np.loadtxt(p['beam_%s_T'%f],unpack=True)
        l,fl_Pol=np.loadtxt(p['beam_%s_pol'%f],unpack=True)

        alm_beamed[0]=hp.sphtfunc.almxfl(alms[0], fl_T)
        alm_beamed[1]=hp.sphtfunc.almxfl(alms[1], fl_Pol)
        alm_beamed[2]=hp.sphtfunc.almxfl(alms[2], fl_Pol)
        
        if pixel=='healpix':
            if pixWin==True:
                pw_T,pw_Pol=hp.sphtfunc.pixwin(nside, pol=True)
                alm_beamed[0]=hp.sphtfunc.almxfl(alm_beamed[0], pw_T)
                alm_beamed[1]=hp.sphtfunc.almxfl(alm_beamed[1], pw_Pol)
                alm_beamed[2]=hp.sphtfunc.almxfl(alm_beamed[2], pw_Pol)

            m=hp.sphtfunc.alm2map(alm_beamed, nside)

        if pixel=='car':
            template= enmap.read_map('%s'%p['template'])
            m=curvedsky.alm2map(alm_beamed, template)
            if pixWin==True:
                m= enmap.apply_window(m, pow=1.0)

        for s in range(nSplit):
            
            maps=iso_map_utils.add_noise(m,f,s,noise,nSplit,pixel)
            if pixel=='car' and pixWin==True:
                maps= enmap.apply_window(maps, pow=-1.0)

            
            fName='split_%d_%s.fits'%(s,f)
            
            iso_map_utils.write_map(maps,mapDir,fName,pixel)


