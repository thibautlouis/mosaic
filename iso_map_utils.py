#!/usr/bin/env python

import healpy as hp
import numpy as np
import pylab as plt
import os
import matplotlib.style
import matplotlib as mpl
mpl.style.use('classic')

def create_directory(dirName):
    try:
        os.makedirs(dirName)
    except:
        pass

def add_noise(m,rmsT,rmsP,nsplit,pixel):
    
    rad_to_arcmin=60*180/np.pi
    m_noisy=m.copy()
    
    if pixel=='healpix':
        nside=hp.pixelfunc.get_nside(m)
        pixArea= hp.pixelfunc.nside2pixarea(nside)*rad_to_arcmin**2
        rmsT=rmsT/np.sqrt(pixArea)
        rmsP=rmsP/np.sqrt(pixArea)

        size=len(m[0])
        m_noisy[0]=m[0]+np.random.randn(size)*np.sqrt(nsplit)*rmsT
        m_noisy[1]=m[1]+np.random.randn(size)*np.sqrt(nsplit)*rmsP
        m_noisy[2]=m[2]+np.random.randn(size)*np.sqrt(nsplit)*rmsP
    if pixel=='car':
        pixArea= m.pixsizemap()*rad_to_arcmin**2
        rmsT=rmsT/np.sqrt(pixArea)
        rmsP=rmsP/np.sqrt(pixArea)
        
        size=m[0].shape
        m_noisy[0]=m[0]+np.random.randn(size[0],size[1])*np.sqrt(nsplit)*rmsT
        m_noisy[1]=m[1]+np.random.randn(size[0],size[1])*np.sqrt(nsplit)*rmsP
        m_noisy[2]=m[2]+np.random.randn(size[0],size[1])*np.sqrt(nsplit)*rmsP

    return(m_noisy)




def write_map(maps,dir,fName,pixel):
    if pixel=='healpix':
        hp.fitsfunc.write_map(dir+fName, maps,overwrite=True)
    if pixel=='car':
        from enlib import enmap
        enmap.write_map(dir+fName, maps)

def read_map(dir,fName,pixel):
    if pixel=='healpix':
        try:
            m=hp.fitsfunc.read_map(dir+fName,field=(0,1,2),verbose=False)
        except IndexError:
            m=hp.fitsfunc.read_map(dir+fName,field=0,verbose=False)
    if pixel=='car':
        from enlib import enmap
        m=enmap.read_map(dir+fName)
    return m

def cut_patch_car(map,ra0,ra1,dec0,dec1):
    box= np.array( [[ dec0, ra1], [dec1, ra0]])*np.pi/180
    map=map.submap(box)
    return(map)
