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

def write_in_hdf5(name,lb,cb_dict):
    spec=file.create_group(name)
    array=iso_ps_utils.get_cl_array(lb,cb_dict)
    spec.create_dataset(name='data',data=array,dtype='float')


def get_spectra(mapDir,auxDir,mcmDir,specDir,winList,nSplits,niter,lmax,binningFile,type,hdf5,pixel,survey_mask_coordinates=None,removeMean=None,thetaCut=None):

    if hdf5==True:
        file = h5py.File('%s.hdf5'%(specDir), 'w')
    else:
        file= None
        iso_map_utils.create_directory(specDir)

    if pixel=='car' and survey_mask_coordinates is not None:
        ra0,ra1,dec0,dec1=np.loadtxt(survey_mask_coordinates,unpack=True, usecols=range(1,5),ndmin=2)

    nPatch,freq=iso_window_utils.get_frequency_list(winList)

    if thetaCut==True:
        theta_range=np.loadtxt('%s/window_theta_range.txt'%(auxDir))

    for i in range(nPatch):
        spec_list = open('%s/spectra_list_%03d.txt'%(auxDir,i),mode="w")
        spec_list_combined = open('%s/spectra_list_combined_%03d.txt'%(auxDir,i),mode="w")

        if len(freq[i]) !=0:
            patch_name='patch_%03d'%i

            alms={}
            for f1 in freq[i]:
                winName_T= 'window_%s_%03d_T.fits'%(f1,i)
                winName_pol= 'window_%s_%03d_pol.fits'%(f1,i)
                window_T=  iso_map_utils.read_map(auxDir,winName_T,pixel)
                window_pol=  iso_map_utils.read_map(auxDir,winName_pol,pixel)
                window=[window_T,window_pol,window_pol]
                del window_T,window_pol
                for s in range(nSplits):
                    fName='split_%d_%s.fits'%(s,f1)
                    map=iso_map_utils.read_map(mapDir,fName,pixel)
                    
                    if pixel=='car' and survey_mask_coordinates is not None:
                        map=iso_map_utils.cut_patch_car(map,ra0[i],ra1[i],dec0[i],dec1[i])
                    
                    if removeMean==True:
                        for k in range(3):
                            map[k]-=np.mean(map[k]*window[k])

                    for ii in range(3):
                        map[ii]*=window[ii]
                    
                    thetas=None
                    if pixel=='healpix' and nPatch!=1 and thetaCut==True:
                        thetas=theta_range[i,:]

                    alms[f1,s]=iso_ps_utils.map2alm(map,pixel,niter,lmax=lmax,theta_range=thetas)
        
            count1=0
            for f1 in freq[i]:
                count2=0
                for f2 in freq[i]:
                    if count2>count1:
                        continue
                    mcm=np.loadtxt('%s/mode_coupling_%03d_%sGHzx%sGHz.dat'%(mcmDir,i,f1,f2))

                    dict_list_auto=[]
                    dict_list_cross=[]

                    count_s1=0
                    for s1 in range(nSplits):
                        count_s2=0
                        for s2 in range(nSplits):
                
                            if ((f1==f2) & (count_s2>count_s1)):
                                break
                        
                            l,cl_dict=iso_ps_utils.get_cl_dict(alms[f1,s1],alms[f2,s2])
                            lb,cb_dict=iso_ps_utils.bin_cl_dict(l,cl_dict,binningFile,lmax,type=type)

                            cb_dict=iso_ps_utils.apply_mcm(cb_dict,mcm,direct_invert=True)
                        
                            fName='%s_%sGHzx%sGHz_split%dxsplit%d_binned'%(type,f1,f2,s1,s2)
                            iso_ps_utils.write_cl_dict(specDir,patch_name,fName,lb,cb_dict,hdf5,file)
                            
                            spec_list.write('%s \n'%fName)

                            
                            if s1 == s2  and f1 == f2:
                                print 'auto',s1,s2,f1,f2
                                dict_list_auto+=[cb_dict]
                            else:
                                print 'cross',s1,s2,f1,f2
                                dict_list_cross+=[cb_dict]

                            count_s2+=1
                        count_s1+=1
                    
                    print '%d auto spectra generated'%len(dict_list_auto)
                    print '%d cross spectra generated'%len(dict_list_cross)
                    
                    cb_cross_mean=iso_ps_utils.get_mean_dict(dict_list_cross)
                    cb_auto_mean=iso_ps_utils.get_mean_dict(dict_list_auto)
                    nb_mean=iso_ps_utils.get_noise_dict(cb_auto_mean,cb_cross_mean)
                    
                    fName='%s_%sGHzx%sGHz_mean_auto_binned'%(type,f1,f2)
                    iso_ps_utils.write_cl_dict(specDir,patch_name,fName,lb,cb_auto_mean,hdf5,file)
                    fName='%s_%sGHzx%sGHz_mean_cross_binned'%(type,f1,f2)
                    iso_ps_utils.write_cl_dict(specDir,patch_name,fName,lb,cb_cross_mean,hdf5,file)
                    spec_list_combined.write('%s \n'%fName)
                    fName='%s_%sGHzx%sGHz_mean_noise_binned'%(type,f1,f2)
                    iso_ps_utils.write_cl_dict(specDir,patch_name,fName,lb,nb_mean,hdf5,file)
                    

                    count2+=1
                count1+=1
                
                spec_list.close()
                spec_list_combined.close()
                    
    if hdf5==True:
        file.close()

def get_spectra_namaster(p,mapDir,auxDir,mcmDir,specDir,winList,nSplits,lmax,type,pixel,nlb=None):
    import pymaster as nmt

    iso_map_utils.create_directory(specDir)
    
    if pixel=='car':
        print '*'*50
        print 'error: naMaster doesnt support car pixelization yet'
        print '*'*50
        sys.exit()
    if type=='Dl':
        print '*'*50
        print 'error: naMaster doesnt Dl yet only Cl'
        print '*'*50
        sys.exit()

    def compute_master(f_a,f_b,wsp) :
        cl_coupled=nmt.compute_coupled_cell(f_a,f_b,n_iter=0)
        cl_decoupled=wsp.decouple_cell(cl_coupled)
        return cl_decoupled

    nPatch,freq=iso_window_utils.get_frequency_list(winList)
    
    for i in range(nPatch):
        spec_list = open('%s/spectra_list_%03d.txt'%(specDir,i),mode="w")

        if len(freq[i]) !=0:
            patch_name='patch_%03d'%i
            
            field_0={}
            field_1={}
            field_2={}

            for f1 in freq[i]:
                winName_T= 'window_%s_%03d_T.fits'%(f1,i)
                winName_pol= 'window_%s_%03d_pol.fits'%(f1,i)
                window_T=  iso_map_utils.read_map(auxDir,winName_T,pixel)
                window_pol=  iso_map_utils.read_map(auxDir,winName_pol,pixel)
                
                nside=hp.pixelfunc.get_nside(window_T)
                ell,beam_T=np.loadtxt(p['beam_%s_T'%f1],unpack=True)
                ell,beam_pol=np.loadtxt(p['beam_%s_Pol'%f1],unpack=True)

                for s in range(nSplits):
                    fName='split_%d_%s.fits'%(s,f1)
                    map=iso_map_utils.read_map(mapDir,fName,pixel)
                    field_0[f1,s]=nmt.NmtField(window_T,[map[0]], beam=beam_T[:3*nside])
                    field_2[f1,s]=nmt.NmtField(window_pol,[map[1],map[2]], beam=beam_pol[:3*nside])
                        
            b=nmt.NmtBin(nside,nlb=nlb)
            lb=b.get_effective_ells()

            count1=0
            for f1 in freq[i]:
                count2=0
                for f2 in freq[i]:
                    if count2>count1:
                        continue
                    
                    w0=nmt.NmtWorkspace()
                    w0.compute_coupling_matrix(field_0[f1,0],field_0[f2,0],b)
                    w1=nmt.NmtWorkspace()
                    w1.compute_coupling_matrix(field_0[f1,0],field_2[f2,0],b)
                    w2=nmt.NmtWorkspace()
                    w2.compute_coupling_matrix(field_2[f1,0],field_0[f2,0],b)
                    w3=nmt.NmtWorkspace()
                    w3.compute_coupling_matrix(field_2[f1,0],field_2[f2,0],b)
                    
                    count_s1=0
                    for s1 in range(nSplits):
                        count_s2=0
                        for s2 in range(nSplits):
                                
                            if ((f1==f2) & (count_s2>count_s1)):
                                break

                            cb_dict={}
                            cb_dict['TT']=compute_master(field_0[f1,s1],field_0[f2,s2],w0)[0]
                            spin1=compute_master(field_0[f1,s1],field_2[f2,s2],w1)
                            cb_dict['TE']=spin1[0]
                            cb_dict['TB']=spin1[1]
                            spin1=compute_master(field_2[f1,s1],field_0[f2,s2],w2)
                            cb_dict['ET']=spin1[0]
                            cb_dict['BT']=spin1[1]
                            spin2=compute_master(field_2[f1,s1],field_2[f2,s2],w3)
                            cb_dict['EE']=spin2[0]
                            cb_dict['EB']=spin2[1]
                            cb_dict['BE']=spin2[2]
                            cb_dict['BB']=spin2[3]

                            fName='%s_%sGHzx%sGHz_split%dxsplit%d_binned'%(type,f1,f2,s1,s2)
                            spec_list.write('%s \n'%fName)

                            iso_ps_utils.write_cl_dict(specDir,patch_name,fName,lb,cb_dict)
                            count_s2+=1
                        count_s1+=1
                    count2+=1
                count1+=1
                
                spec_list.close()

