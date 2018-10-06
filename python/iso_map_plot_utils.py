#!/usr/bin/env python

import healpy as hp
import numpy as np
import pylab as plt
import iso_map_utils
import iso_window_utils
import iso_ps_utils
import os
import sys
import h5py


#This plotting routine need serious rewriting, it's pretty ugly

def theta_phi_healpix(nside):
    n_pts=12*nside**2
    theta=np.zeros(n_pts)
    phi=np.zeros(n_pts)
    for i in range(n_pts):
        theta[i],phi[i] = hp.pixelfunc.pix2ang(nside, i)
    return(theta,phi)

def lon_lat_healpix(nside):
    theta,phi=theta_phi_healpix(nside)
    lat=(np.pi/2-theta)*180/np.pi
    long=phi*180/np.pi
    return(long,lat)

def get_disk(nside,radius,theta,phi):
    pix=hp.get_interp_weights(nside,theta,phi)[0][0]
    vec=hp.pix2vec(nside,pix)
    disk=hp.query_disc(nside, vec, radius)
    return(disk)

def plot_maps(maps,pixel,mask=0,color='planck',color_range=None,png_file=None,gnomview_coord=None):
    
    def get_xsize_from_radius(radius):
        reso=1.5
        radius*=60/reso
        xsize=2*radius
        return(xsize)
    
    if pixel=='healpix':
        
        if color=='planck':
            from matplotlib.colors import ListedColormap
            path=(os.path.dirname(os.path.realpath(__file__)))
            colombi1_cmap = ListedColormap(np.loadtxt('%s/../data/Planck_Parchment_RGB.txt'%path)/255.)
            colombi1_cmap.set_bad("white") # color of missing pixels
            colombi1_cmap.set_under("white")
            cmap = colombi1_cmap
        else:
            cmap= plt.get_cmap('jet')

        if maps.ndim==1:
            
            if gnomview_coord is None:
                if color_range is not None:
                    hp.mollview(maps,min=color_range,max=color_range,cmap=cmap, notext=True)
                else:
                    hp.mollview(maps,cmap=cmap, notext=True)

                if png_file is not None:
                    plt.savefig(png_file+'.png')
                    plt.clf()
                    plt.close
                else:
                    plt.show()
    
            else:
                xsize=get_xsize_from_radius(gnomview_coord[2])
                if color_range is not None:
                    hp.visufunc.gnomview(maps,min=-color_range,max=color_range, rot=(gnomview_coord[0],gnomview_coord[1]),xsize=xsize,cmap=cmap)
                else:
                    hp.visufunc.gnomview(maps, rot=(gnomview_coord[0],gnomview_coord[1]),xsize=xsize,cmap=cmap)
                if png_file is not None:
                    plt.savefig(png_file+'_projected.png')
                    plt.clf()
                    plt.close
                else:
                    plt.show()
        else:
            fields=['T','Q','U']
            if color_range is not None:
                min={}
                max={}
                for i,l1 in enumerate(fields):
                    min[l1]=-color_range[i]
                    max[l1]=color_range[i]
            
            for map,l1 in zip(maps,fields):
                
                if gnomview_coord is None:
                    if color_range is not None:
                        hp.mollview(map,min=min[l1],max=max[l1],cmap=cmap, notext=True)
                    else:
                        hp.mollview(map,cmap=cmap, notext=True)

                    if png_file is not None:
                        plt.savefig(png_file+'_%s'%l1+'.png')
                        plt.clf()
                        plt.close
                    else:
                        plt.show()
                else:
                    xsize=get_xsize_from_radius(gnomview_coord[2])
                    if color_range is not None:
                        hp.visufunc.gnomview(map,min=min[l1],max=max[l1], rot=(gnomview_coord[0],gnomview_coord[1]),xsize=xsize,cmap=cmap)
                    else:
                        hp.visufunc.gnomview(map, rot=(gnomview_coord[0],gnomview_coord[1]),xsize=xsize,cmap=cmap)

                    if png_file is not None:
                        plt.savefig(png_file+'_%s'%l1+'_projected.png')
                        plt.clf()
                        plt.close
                    else:
                        plt.show()
    if pixel=='car':
        from enlib import enplot
    
        if maps.ndim==2:
            plots = enplot.get_plots(maps,mask=mask,color=color)
            for plot in plots:
                if png_file is not None:
                    enplot.write_plot(png_file+'.png', plot)
                else:
                    plot.img.show()
        if maps.ndim==3:
            fields=['T','Q','U']

            if color_range is not None:
                colors='%s:%s:%s'%(color_range[0],color_range[1],color_range[2])
                plots = enplot.get_plots(maps,mask=mask,color=color,range=colors,colorbar=1)
            else:
                plots = enplot.get_plots(maps,mask=mask,color=color)

            for (plot,l1) in zip(plots,fields):
                if png_file is not None:
                    enplot.write_plot(png_file+'_%s'%l1+'.png', plot)
                else:
                    plot.img.show()

def plot_all_windows(auxDir,plotDir,pixel,winList,tessel_healpix=None):
    
    if pixel=='healpix':
        if tessel_healpix['apply']:
            long_c,lat_c=lon_lat_healpix(tessel_healpix['patch_nside'])
            radius=tessel_healpix['radius']


    nPatch,freq=iso_window_utils.get_frequency_list(winList)
    for i in range(nPatch):
        for f1 in freq[i]:
            winName_T= 'window_%s_%03d_T'%(f1,i)
            winName_pol= 'window_%s_%03d_pol'%(f1,i)
            window_T=  iso_map_utils.read_map(auxDir,winName_T+'.fits',pixel)
            window_pol=  iso_map_utils.read_map(auxDir,winName_pol+'.fits',pixel)
            
            gnomview_coord=None
            if pixel=='healpix':
                if tessel_healpix['apply']:
                    gnomview_coord=[long_c[i],lat_c[i],radius]

            plot_maps(window_T,pixel, png_file=plotDir+'/'+winName_T,gnomview_coord=gnomview_coord)
            plot_maps(window_pol,pixel, png_file=plotDir+'/'+winName_pol,gnomview_coord=gnomview_coord)

def plot_all_maps(auxDir,mapDir,plotDir,pixel,winList,nSplits,color_range=None,survey_mask_coordinates=None,tessel_healpix=None):
    
    if pixel=='car':
        ra0,ra1,dec0,dec1=np.loadtxt(survey_mask_coordinates,unpack=True, usecols=range(1,5),ndmin=2)
    
    if pixel=='healpix':
        if tessel_healpix['apply']:
            long_c,lat_c=lon_lat_healpix(tessel_healpix['patch_nside'])
            radius=tessel_healpix['radius']
        else:
            gnomview_coord=None

    nPatch,freq=iso_window_utils.get_frequency_list(winList)

    for i in range(nPatch):
        for f1 in freq[i]:
            winName_T= 'window_%s_%03d_T'%(f1,i)
            winName_pol= 'window_%s_%03d_pol'%(f1,i)
            window_T=  iso_map_utils.read_map(auxDir,winName_T+'.fits',pixel)
            window_pol=  iso_map_utils.read_map(auxDir,winName_pol+'.fits',pixel)
            window=[window_T,window_pol,window_pol]
            del window_T,window_pol
            
            for s in range(nSplits):
                fName='split_%d_%s'%(s,f1)
                maps=iso_map_utils.read_map(mapDir,fName+'.fits',pixel)
                if pixel=='car':
                    maps=iso_map_utils.cut_patch_car(maps,ra0[i],ra1[i],dec0[i],dec1[i])
                for ii in range(3):
                    maps[ii]*=window[ii]
                gnomview_coord=None
                if pixel=='healpix':

                    if tessel_healpix['apply']:
                        gnomview_coord=[long_c[i],lat_c[i],radius]
                plot_maps(maps,pixel,mask=0,color='planck',color_range=color_range,png_file=plotDir+'/'+fName+'_%03d'%i,gnomview_coord=gnomview_coord)

def plot_survey_map(p,auxDir,mapDir,plotDir,pixel,winList,color_range=None,survey_mask_coordinates=None,tessel_healpix=None):
    if pixel=='car':
        from enlib import enmap,enplot
        ra0,ra1,dec0,dec1=np.loadtxt(survey_mask_coordinates,unpack=True, usecols=range(1,5),ndmin=2)
        nPatch,freq=iso_window_utils.get_frequency_list(winList)
        for i in range(nPatch):
            for f1 in freq[i]:
                fName='split_%d_%s'%(0,f1)
                maps=iso_map_utils.read_map(mapDir,fName+'.fits',pixel)
                T = enmap.downgrade(maps[0], 40)
                img = enplot.draw_map_field(T, enplot.parse_args("a -r 1000 -g -m 0"))
                img_box = enplot.draw_annotations(T, [["rect", dec0[i], ra0[i], 0, 0, dec1[i], ra1[i], 0, 0]], None)
                img_tot = enplot.merge_images([img, img_box])
                img_tot.save('%s/im_%03d.png'%(plotDir,i))
    if pixel=='healpix':
        mask=p['mask']
        deg_to_rad=np.pi/180
        if tessel_healpix['apply']:
            theta_c,phi_c=theta_phi_healpix(tessel_healpix['patch_nside'])
            nPatch,freq=iso_window_utils.get_frequency_list(winList)
            for i in range(nPatch):
                for f1 in freq[i]:
                    fName='split_%d_%s'%(0,f1)
                    maps=iso_map_utils.read_map(mapDir,fName+'.fits',pixel)
                    nside=hp.pixelfunc.get_nside(maps)
                    T=maps[0]
                    
                    if mask['gal_mask_T']==True:
                        mName=p['gal_mask_T_%s_%s'%(pixel,f1)]
                        T=iso_window_utils.apply_mask(T,mName,pixel,box=None)

                    template=T.copy()
                    template*=0
                    radius1=tessel_healpix['radius']*deg_to_rad
                    disk1=get_disk(nside,radius1,theta_c[i],phi_c[i])
                    template[disk1]+=1
                    disk2=get_disk(nside,radius1*1.05,theta_c[i],phi_c[i])
                    template[disk2]+=1
                    id=np.where(template==1)
                    T[id]=np.nan
                    plot_maps(T,pixel,png_file=plotDir+'/im_%03d_%s.png'%(i,f1))
