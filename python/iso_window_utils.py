#!/usr/bin/env python

import healpy as hp
import numpy as np
import pylab as plt
import iso_map_utils
import iso_apodization_utils
import os

def get_frequency_list(win_list):
    winList=np.genfromtxt(win_list,unpack=True,dtype='str')
    num = [int(i) for i in winList[0]]
    nPatch=np.max(num)+1
    if nPatch==1:
        return(nPatch,[winList[1:]])
    freq={}
    for i in range(nPatch):
        freq[i]= list(np.transpose(winList)[i][1:])
        try:
            for j in range(len(freq[i])):
                freq[i].remove('000')
        except:
            pass
    
    return(nPatch,freq)

def get_spectra_list(win_list):
    
    nPatch,freq=get_frequency_list(win_list)
    spec_list=[]
    for i in range(nPatch):
        if len(freq[i]) !=0:
            count1=0
            for f1 in freq[i]:
                count2=0
                for f2 in freq[i]:
                    if count2>count1:
                        continue
                    spec_list+=[[i,'%s'%f1,'%s'%f2]]
                    count2+=1
                count1+=1
    return spec_list

def apply_mask(map,mName,pixel,box):
    
    if pixel=='healpix':
        nside=hp.pixelfunc.get_nside(map)
        m,header=hp.fitsfunc.read_map(mName,h=True,verbose=False)
        mask_nside=header[13][1]
        if mask_nside != nside:
            raise ValueError, "mask nside and map nside disagree"
    if pixel=='car':
        from enlib import enmap
        m=enmap.read_map(mName)
        if box is not None:
            m=m.submap(box)
    map*=m
    return(map)

def apply_ps_mask(mask,template,f,pixel,box=None):
    
    maskT=template.copy()
    maskPol=template.copy()
    if mask['point_mask_T']==True:
        mName=mask['pts_mask_T_%s_%s'%(pixel,f)]
        maskT=apply_mask(maskT,mName,pixel,box)
    if mask['point_mask_pol']==True:
        mName=mask['pts_mask_pol_%s_%s'%(pixel,f)]
        maskPol=apply_mask(maskPol,mName,pixel,box)
        
    return(maskT,maskPol)

def apply_gal_mask(mask,template,f,pixel,box=None):
    
    maskT=template.copy()
    maskPol=template.copy()
    if mask['gal_mask_T']==True:
        mName=mask['gal_mask_T_%s_%s'%(pixel,f)]
        maskT=apply_mask(maskT,mName,pixel,box)
    if mask['gal_mask_pol']==True:
        mName=mask['gal_mask_pol_%s_%s'%(pixel,f)]
        maskPol=apply_mask(maskPol,mName,pixel,box)
    
    return(maskT,maskPol)

def apply_missing_pixel_mask(mask,template,f,pixel,box=None):
    
    maskT=template.copy()
    maskPol=template.copy()
    if mask['missing_pixel_mask_T']==True:
        mName=mask['missing_pixel_mask_T_%s_%s'%(pixel,f)]
        maskT=apply_mask(maskT,mName,pixel,box)
    if mask['missing_pixel_mask_pol']==True:
        mName=mask['missing_pixel_mask_pol_%s_%s'%(pixel,f)]
        maskPol=apply_mask(maskPol,mName,pixel,box)
    
    return(maskT,maskPol)



def mask_healpix(auxDir, freqTags, pixel, tessel_healpix, mask, apo_type, apo_radius, template, plot=False):
    
    deg_to_rad=np.pi/180

    def get_fsky(mask):
        return(np.sum(mask)/np.sum(mask*0+1))

    survey_mask=template*0+1
    
    ps_maskT={}
    ps_maskPol={}
    missing_pixel_maskT={}
    missing_pixel_maskPol={}
    for f in freqTags:
        ps_maskT[f],ps_maskPol[f]=apply_ps_mask(mask,survey_mask,f,pixel)
        missing_pixel_maskT[f],missing_pixel_maskPol[f]=apply_missing_pixel_mask(mask,survey_mask,f,pixel)

        ps_maskT[f]=iso_apodization_utils.create_apodization(ps_maskT[f], pixel, apo_type,apo_radius['ps'])
        ps_maskPol[f]=iso_apodization_utils.create_apodization(ps_maskPol[f], pixel, apo_type,apo_radius['ps'])
        missing_pixel_maskT[f]=iso_apodization_utils.create_apodization(missing_pixel_maskT[f], pixel, apo_type,apo_radius['miss_pix'])
        missing_pixel_maskPol[f]=iso_apodization_utils.create_apodization(missing_pixel_maskPol[f], pixel, apo_type,apo_radius['miss_pix'])

    win_list = open('%s/window_list.txt'%auxDir,mode="w")

    if tessel_healpix['apply']:
        
        def theta_phi_healpix(nside):
            n_pts=12*nside**2
            theta=np.zeros(n_pts)
            phi=np.zeros(n_pts)
            for i in range(n_pts):
                theta[i],phi[i] = hp.pixelfunc.pix2ang(nside, i)
            return(theta,phi)

        def get_disk(nside,radius,theta,phi):
            pix=hp.get_interp_weights(nside,theta,phi)[0][0]
            vec=hp.pix2vec(nside,pix)
            disk=hp.query_disc(nside, vec, radius)
            return(disk)

        theta_c,phi_c=theta_phi_healpix(tessel_healpix['patch_nside'])
        
        theta_min=theta_c-(tessel_healpix['radius']+1)*deg_to_rad
        theta_max=theta_c+(tessel_healpix['radius']+1)*deg_to_rad
        
        np.savetxt('%s/window_theta_range.txt'%auxDir,np.transpose([theta_min,theta_max]))
                   
        n_disks=len(theta_c)
        nside=hp.pixelfunc.get_nside(template)

        print 'generate %d patch from the Healpix map'%n_disks
        
        for i in range(n_disks):
            radius=tessel_healpix['radius']*deg_to_rad
            survey_mask=template*0
            disk=get_disk(nside,radius,theta_c[i],phi_c[i])
            survey_mask[disk]=1
            fsky_disk=get_fsky(survey_mask)

            win_list.write('%03d'%i)

            for f in freqTags:
                surveymask_T,surveymask_Pol=apply_gal_mask(mask,survey_mask,f,pixel)
                fsky_mask=get_fsky(surveymask_T)
                use_frac=fsky_mask/fsky_disk
                
                if fsky_mask/fsky_disk >tessel_healpix['cut_threshold'] :
                    #patch_list.write('%03d %0.4f %0.4f %s %0.4f \n'%(i,phi_c[i],theta_c[i],f,use_frac))
                    
                    surveymask_T=iso_apodization_utils.create_apodization(surveymask_T, pixel, apo_type,apo_radius['survey'])
                    surveymask_Pol=iso_apodization_utils.create_apodization(surveymask_Pol, pixel, apo_type,apo_radius['survey'])

                    surveymask_T*=ps_maskT[f]
                    surveymask_Pol*=ps_maskPol[f]
                    surveymask_T*=missing_pixel_maskT[f]
                    surveymask_Pol*=missing_pixel_maskPol[f]

                    if plot:
                        iso_map_utils.plot(surveymask_T,pixel)
                        iso_map_utils.plot(surveymask_Pol,pixel)
                
                    win_list.write(' %s '%f)

                    winName='window_%s_%03d'%(f,i)

                    iso_map_utils.write_map(surveymask_T,auxDir,winName+'_T.fits',pixel)
                    iso_map_utils.write_map(surveymask_Pol,auxDir,winName+'_pol.fits',pixel)


                else:
                    win_list.write(' 000 ')
                    continue

            win_list.write('\n')

    else:
        survey_mask=template*0+1
        win_list.write('000')

        for f in freqTags:
            surveymask_T,surveymask_Pol=apply_gal_mask(mask,survey_mask,f,pixel)
            
            surveymask_T=iso_apodization_utils.create_apodization(surveymask_T, pixel, apo_type,apo_radius['survey'])
            surveymask_Pol=iso_apodization_utils.create_apodization(surveymask_Pol, pixel, apo_type,apo_radius['survey'])

            surveymask_T*=ps_maskT[f]
            surveymask_Pol*=ps_maskPol[f]

            if plot:
                iso_map_utils.plot(surveymask_T,pixel)
                iso_map_utils.plot(surveymask_Pol,pixel)
            
            winName='window_%s_000'%(f)
            win_list.write(' %s '%f)

            iso_map_utils.write_map(surveymask_T,auxDir,winName+'_T.fits',pixel)
            iso_map_utils.write_map(surveymask_Pol,auxDir,winName+'_pol.fits',pixel)


    win_list.close()

def rectangle_mask_car(auxDir, freqTags, pixel, survey_mask_coordinate, mask, apo_type, apo_radius, template, plot=False):
    from enlib import enmap
    
    ra0,ra1,dec0,dec1=np.loadtxt(survey_mask_coordinate,unpack=True, usecols=range(1,5),ndmin=2)
    nPatch=len(ra0)
    
    win_list = open('%s/window_list.txt'%auxDir,mode="w")

    for i in range(nPatch):
        survey_mask=template*0+1

        box= np.array( [[ dec0[i], ra1[i]], [dec1[i], ra0[i] ]])*np.pi/180
        survey_mask=iso_map_utils.cut_patch_car(survey_mask,ra0[i],ra1[i],dec0[i],dec1[i])

        win_list.write('%03d'%i)

        for f in freqTags:
            survey_mask[:,:]=1
            
            ps_maskT,ps_maskPol=apply_ps_mask(mask,survey_mask,f,pixel,box)
            ps_maskT=iso_apodization_utils.create_apodization(ps_maskT, pixel, apo_type,apo_radius['ps'])
            ps_maskPol=iso_apodization_utils.create_apodization(ps_maskPol, pixel, apo_type,apo_radius['ps'])

            survey_mask[:,:]=0
            survey_mask[1:-1,1:-1]=1

            surveymask_T,surveymask_Pol=apply_gal_mask(mask,survey_mask,f,pixel,box)
            surveymask_T=iso_apodization_utils.create_apodization(surveymask_T, pixel, apo_type,apo_radius['survey'])
            surveymask_Pol=iso_apodization_utils.create_apodization(surveymask_Pol, pixel, apo_type,apo_radius['survey'])
            
            surveymask_T*=ps_maskT
            surveymask_Pol*=ps_maskPol

            if plot:
                iso_map_utils.plot(surveymask_T,pixel,mask=None)
                iso_map_utils.plot(surveymask_Pol,pixel,mask=None)


            winName='window_%s_%03d'%(f,i)
            iso_map_utils.write_map(surveymask_T,auxDir,winName+'_T.fits',pixel)
            iso_map_utils.write_map(surveymask_Pol,auxDir,winName+'_pol.fits',pixel)
            win_list.write(' %s '%f)
                
        win_list.write('\n')


    win_list.close()

