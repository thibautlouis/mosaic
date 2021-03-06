#!/usr/bin/env python

import healpy as hp
import numpy as np
import pylab as plt
import iso_map_utils
import os
import scipy.ndimage
from pixell import enmap

def get_distance(binary,pixel):
    
    def write_dict_file():
        file = open("distance.dict",'w')
        file.write("mask_file=tempmask \n")
        file.write("hole_min_size=0 \n")
        file.write("hole_min_surf_arcmin2=0.0 \n")
        file.write("filled_file='' \n")
        file.write("distance_file=tempfile \n")
        file.close()
        return
    
    if pixel=='healpix':
        try:
            os.system('rm -rf tempfile')
            os.system('rm -rf tempmask')
        except:
            pass
        
        hp.fitsfunc.write_map('tempmask', binary)
        write_dict_file()
        os.system('my_process_mask')
        dist=hp.fitsfunc.read_map('tempfile')
        os.system('rm -rf tempfile')
        os.system('rm -rf tempmask')

        dist*=180/np.pi

    if pixel=='car':
        pixSize_arcmin= np.sqrt(binary.pixsize()*(60*180/np.pi)**2)
        print 'pixSize= %0.2f'%pixSize_arcmin
        dist= scipy.ndimage.distance_transform_edt(binary)
        dist*=pixSize_arcmin/60
        dist= enmap.samewcs(dist, binary)
    
    return dist

def apod_C1(binary,radius,pixel):
    
    if radius==0:
        return binary
    else:
        dist=get_distance(binary,pixel)

        if np.max(dist)==0:
            return (dist*0+1)

        win=dist.copy()
        id=np.where(dist> radius)
        win=dist/radius-np.sin(2*np.pi*dist/radius)/(2*np.pi)
        win[id]=1
    
        return(win)

def apod_C2(binary,radius,pixel):
    
    if radius==0:
        return binary
    else:

        dist=get_distance(binary,pixel)
    
        if np.max(dist)==0:
            return (dist*0+1)

        win=dist.copy()
        id=np.where(dist> radius)
        win=1./2-1./2*np.cos(-np.pi*dist/radius)
        win[id]=1
        return(win)


def apod_C3(binary,radius,pixel):
    
    if radius==0:
        return binary
    else:
        if pixel=='healpix':
            return 'C3 not available for healpix pixellisation'
            sys.exit()

        print radius
        win=car_cosine_window(binary,radius,0)
        return(win)

def car_cosine_window(binary,lenApod,pad):
    
    shape= binary.shape
    wcs= binary.wcs
    Ny,Nx=shape
    pixScaleY,pixScaleX= enmap.pixshape(shape, wcs)
    win=binary.copy()
    win=win*0+1
    winX=win.copy()
    winY=win.copy()
    Id=np.ones((Ny,Nx))
    
    degToPix_x=np.pi/180/pixScaleX
    degToPix_y=np.pi/180/pixScaleY
    
    lenApod_x=int(lenApod*degToPix_x)
    lenApod_y=int(lenApod*degToPix_y)
    pad_x=int(pad*degToPix_x)
    pad_y=int(pad*degToPix_y)
    
    for i in range(pad_x,lenApod_x+pad_x):
        r=float(i)-pad_x
        winX[:,i]=1./2*(Id[:,i]-np.cos(-np.pi*r/lenApod_x))
        winX[:,Nx-i-1]=winX[:,i]
    for j in range(pad_y,lenApod_y+pad_y):
        r=float(j)-pad_y
        winY[j,:]=1./2*(Id[j,:]-np.cos(-np.pi*r/lenApod_y))
        winY[Ny-j-1,:]=winY[j,:]
    
    win[:]=0
    win[pad_y:Ny-pad_y,pad_x:Nx-pad_x]+=winX[pad_y:Ny-pad_y,pad_x:Nx-pad_x]*winY[pad_y:Ny-pad_y,pad_x:Nx-pad_x]
    return(win)


def create_apodization(binary, pixel, apo_type, apo_radius):
    
    if apo_type=='C1':
        window=apod_C1(binary,apo_radius,pixel)
    if apo_type=='C2':
        window=apod_C2(binary,apo_radius,pixel)
    if apo_type=='C3':
        window=apod_C3(binary,apo_radius,pixel)

    if pixel=='healpix':
        if apo_type=='C1_namaster':
            import pymaster as nmt
            window=nmt.mask_apodization(binary,apo_radius,apotype="C1")
        if apo_type=='C2_namaster':
            import pymaster as nmt
            window=nmt.mask_apodization(binary,apo_radius,apotype="C2")
        if apo_type=='smooth_namaster':
            import pymaster as nmt
            window=nmt.mask_apodization(binary,apo_radius,apotype="Smooth")

#iso_map_utils.plot(window,pixel,mask=None)

    return window


