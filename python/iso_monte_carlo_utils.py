#!/usr/bin/env python

import healpy as hp
import numpy as np
import pylab as plt
import iso_spectra_utils
import iso_map_utils
import iso_window_utils
import iso_ps_utils
import iso_dict
import sys
import h5py
import os

def bintheory(Bbl_T,Bbl_cross,Bbl_pol,clth):
    cbth={}
    cbth['TT']=np.dot(Bbl_T,clth['TT'][:lmax])
    cbth['TE']=np.dot(Bbl_cross,clth['TE'][:lmax])
    cbth['ET']=cbth['TE'][:lmax]
    cbth['TB']=cbth['TT'][:lmax]*0
    cbth['BT']=cbth['TT'][:lmax]*0
    
    vec=np.zeros((3*lmax))
    vec[:lmax]=clth['EE'][:lmax]
    vec[lmax:2*lmax]=clth['EE'][:lmax]*0
    vec[2*lmax:3*lmax]=clth['BB'][:lmax]
    
    Nbin= Bbl_pol.shape[0]/3
    vecb=np.dot(Bbl_pol,vec)
    cbth['EE']=vecb[:Nbin]
    cbth['EB']=vecb[Nbin:2*Nbin]
    cbth['BB']=vecb[2*Nbin:3*Nbin]
    cbth['BE']=cbth['EB']
    return(cbth)

def get_mean_and_std(specDir,mcDir,iStart,iStop,winList,nSplits,type,hdf5):
    
    fields=['TT','EE','BB','TE','EB','TB','ET','BE','BT']

    if hdf5==True:
        file_mc = h5py.File('%s.hdf5'%(mcDir), 'w')
    else:
        file_mc= None
        iso_map_utils.create_directory(mcDir)

    nPatch,freq=iso_window_utils.get_frequency_list(winList)
    
    for i in range(nPatch):
        
        patchName='patch_%03d'%i
        
        if len(freq[i]) !=0:
            count1=0
            for f1 in freq[i]:
                count2=0
                for f2 in freq[i]:
                    if count2>count1:
                        continue
                    count_s1=0
                    for s1 in range(nSplits):
                        count_s2=0
                        for s2 in range(nSplits):
                            if ((f1==f2) & (count_s2>count_s1)):
                                break
                                        
                            cb_dict_list=[]
                            for iii in range(iStart,iStop):
                                simspecDir=specDir+'_%03d'%iii
                                if hdf5==True:
                                    file = h5py.File('%s.hdf5'%(simspecDir), 'r')
                                else:
                                    file=None

                                fName='%s_%sGHzx%sGHz_split%dxsplit%d_binned'%(type,f1,f2,s1,s2)
                                lb,cb_dict=iso_ps_utils.read_cl_dict(simspecDir,patchName,fName,hdf5,file)

                                cb_dict_list+=[cb_dict]

                            mean,std=iso_ps_utils.get_mean_dict(cb_dict_list,get_std=True)
                            meanName='mc_averaged_'+ fName
                            iso_ps_utils.write_cl_dict(mcDir,patchName,meanName,lb,mean,hdf5,file_mc)
                            stdName='std_'+ fName
                            iso_ps_utils.write_cl_dict(mcDir,patchName,stdName,lb,std,hdf5,file_mc)

                                    
                            count_s2+=1
                        count_s1+=1
                    
                    cb_dict_list_1=[]
                    cb_dict_list_2=[]
                    cb_dict_list_3=[]
                        
                    fName_1='%s_%sGHzx%sGHz_mean_auto_binned'%(type,f1,f2)
                    fName_2='%s_%sGHzx%sGHz_mean_cross_binned'%(type,f1,f2)
                    fName_3='%s_%sGHzx%sGHz_mean_noise_binned'%(type,f1,f2)

                        
                    for iii in range(iStart,iStop):
                        simspecDir=specDir+'_%03d'%iii
                        if hdf5==True:
                            file = h5py.File('%s.hdf5'%(simspecDir), 'r')
                        else:
                            file=None
                            
                        lb,cb_dict_1=iso_ps_utils.read_cl_dict(simspecDir,patchName,fName_1,hdf5,file)
                        lb,cb_dict_2=iso_ps_utils.read_cl_dict(simspecDir,patchName,fName_2,hdf5,file)
                        lb,cb_dict_3=iso_ps_utils.read_cl_dict(simspecDir,patchName,fName_3,hdf5,file)

                        cb_dict_list_1+=[cb_dict_1]
                        cb_dict_list_2+=[cb_dict_2]
                        cb_dict_list_3+=[cb_dict_3]

                    mean,std=iso_ps_utils.get_mean_dict(cb_dict_list_1,get_std=True)
                    meanName='mc_averaged_'+ fName_1
                    iso_ps_utils.write_cl_dict(mcDir,patchName,meanName,lb,mean,hdf5,file_mc)
                    stdName='std_'+ fName_1
                    iso_ps_utils.write_cl_dict(mcDir,patchName,stdName,lb,std,hdf5,file_mc)

                    mean,std=iso_ps_utils.get_mean_dict(cb_dict_list_2,get_std=True)
                    meanName='mc_averaged_'+ fName_2
                    iso_ps_utils.write_cl_dict(mcDir,patchName,meanName,lb,mean,hdf5,file_mc)
                    stdName='std_'+ fName_2
                    iso_ps_utils.write_cl_dict(mcDir,patchName,stdName,lb,std,hdf5,file_mc)

                    mean,std=iso_ps_utils.get_mean_dict(cb_dict_list_3,get_std=True)
                    meanName='mc_averaged_'+ fName_3
                    iso_ps_utils.write_cl_dict(mcDir,patchName,meanName,lb,mean,hdf5,file_mc)
                    stdName='std_'+ fName_3
                    iso_ps_utils.write_cl_dict(mcDir,patchName,stdName,lb,std,hdf5,file_mc)

                    count2+=1
                count1+=1
            



