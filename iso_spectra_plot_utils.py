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


def dict_from_theory_file(fields,theoryFile,lmax):
    clth_dict={}
    lth, clth_dict['TT'], clth_dict['EE'], clth_dict['BB'], clth_dict['TE'] = np.loadtxt(theoryFile,unpack=True)
    clth_dict['ET']=clth_dict['TE']
    clth_dict['TB']=clth_dict['TT']*0
    clth_dict['BT']=clth_dict['TB']
    clth_dict['EB']=clth_dict['TT']*0
    clth_dict['BE']=clth_dict['EB']
    for l1 in fields:
        clth_dict[l1]=clth_dict[l1][:lmax]
    lth=lth[:lmax]

    return lth,clth_dict

def get_nlth(noiseInfo,beamName,fields,lmax):
             
    l,fell_T=np.loadtxt(beamName[0],unpack=True)
    l,fell_Pol=np.loadtxt(beamName[1],unpack=True)
    
    rmsT = noiseInfo[1]
    rmsP = noiseInfo[2]
    nl_th={}
    label_th={}
    for l1 in fields:
        nl_th[l1]=np.zeros(lmax)
        label_th[l1]=''
    
    nl_th['TT']=np.ones(lmax)*(rmsT*np.pi/(60*180))**2/fell_T[:lmax]**2*noiseInfo[0]
    label_th['TT']='rms T %0.2f uk.arcmin'%rmsT
    nl_th['EE']=np.ones(lmax)*(rmsP*np.pi/(60*180))**2/fell_Pol[:lmax]**2*noiseInfo[0]
    nl_th['BB']=np.ones(lmax)*(rmsP*np.pi/(60*180))**2/fell_Pol[:lmax]**2*noiseInfo[0]
    label_th['EE']='rms pol %0.2f uk.arcmin'%rmsP
    label_th['BB']='rms pol %0.2f uk.arcmin'%rmsP

    return(nl_th,label_th)

def plot_cl_dict(plotDir,patchName,fName,type,lb,dictList,dictNames,fields,theoryFile,lmax,stdList=None,Bblfile=None,semilog_T=False,semilog_pol=False,noiseInfo=None,beamName=None):
    
    if noiseInfo is not None and beamName is not None:
        nl_th,label_th=get_nlth(noiseInfo,beamName,fields,lmax)

    if theoryFile is not None:
        lth,clth_dict=dict_from_theory_file(fields,theoryFile,lmax)

    if Bblfile is not None:
        Bbl=np.load('%s'%Bblfile)
        cbth_dict=bin_dict(Bbl,plot_cl_dict)

    if type=='Dl':
        fb=lb*0+1
    if type=='Cl':
        fb=lb**2/(2*np.pi)

    count=1
    plt.figure(figsize=(24,12))
    for l1 in fields:
        plt.subplot(3,3,count)
        if l1=='TT' and semilog_T:
            plt.semilogy()
        if l1=='EE' and semilog_pol:
            plt.semilogy()
        if l1=='BB' and semilog_pol:
            plt.semilogy()
        plt.errorbar(lth,clth_dict[l1],label='CMB')
        if noiseInfo is not None:
            plt.errorbar(lth,clth_dict[l1]+nl_th[l1]*lth*(lth+1)/(2*np.pi))
            plt.errorbar(lth,nl_th[l1]*lth*(lth+1)/(2*np.pi))
            plt.legend()
        if Bblfile:
            plt.errorbar(lb,cbth_dict[l1])
        
        if stdList is not None:
            for dict,std,name in zip(dictList,stdList,dictNames):
                plt.errorbar(lb,dict[l1]*fb,std[l1]*fb,fmt='o',label='%s'%name)
        else:
            for dict,name in zip(dictList,dictNames):
                if name=='namaster':
                    plt.errorbar(lb,dict[l1]*fb,fmt='-',label='%s'%name)
                else:
                    plt.errorbar(lb,dict[l1]*fb,fmt='o',label='%s'%name)

        plt.ylabel(r'$D^{%s}_{\ell}$'%l1,fontsize=18)
        if count==1:
            plt.legend()
        if count>6:
            plt.xlabel(r'$\ell$',fontsize=22)
        count+=1

    plt.suptitle('%s %s'%(patchName,fName),fontsize=22)
    plt.savefig('%s/%s_%s.png'%(plotDir,patchName,fName),bbox_inches='tight')
    plt.clf()
    plt.close()


def plot_all_spectra(mcmDir,specDir,plotDir,winList,nSplits,hdf5,type,theoryFile,lmax,compare_mosaic_namaster=False,white_noise_level=None,mcDir=None):


    def increment_dict_list(specDir,fName,name, cbList, dictNamesList ,stdList):

        lb,cb_dict=iso_ps_utils.read_cl_dict(specDir,patchName,fName,hdf5,file)
        cbList+=[cb_dict]
        dictNamesList+=[name]
        if mcDir is not None:
            lb,std_dict=iso_ps_utils.read_cl_dict(mcDir,patchName,'std_'+fName,hdf5,file)
            stdList+=[std_dict]
        else:
            stdList=None
        return (lb, cbList, dictNamesList,stdList)


    fields=['TT','EE','BB','TE','EB','TB','ET','BE','BT']

    if hdf5==True:
        file = h5py.File('%s.hdf5'%(specDir), 'r')
        if compare_mosaic_namaster==True:
            file_nm = h5py.File('%s'%(specDir)+'_namaster.hdf5', 'r')
    else:
        file=None
        file_nm=None

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

                            cbList=[]
                            dictNamesList=[]
                            
                            if mcDir is not None:
                                stdList=[]
                            else:
                                stdList=None
                            
                            fName='%s_%sGHzx%sGHz_split%dxsplit%d_binned'%(type,f1,f2,s1,s2)
                            lb, cbList, dictNamesList,stdList= increment_dict_list(specDir,fName,'mosaic',cbList, dictNamesList ,stdList)

                            if compare_mosaic_namaster:
                                id= specDir[-4:]
                                naMasterSpecDir=specDir[:-4]+'_namaster'+id
                                lb, cbList, dictNamesList,stdList=increment_dict_list(naMasterSpecDir,fName,'namaster',cbList, dictNamesList ,stdList)

                            plot_cl_dict(plotDir,patchName,fName,type,lb,cbList,dictNamesList,fields,theoryFile,lmax,stdList=stdList,Bblfile=None,semilog_T=True,semilog_pol=True)
          
                            
                            count_s2+=1
                        count_s1+=1
                        
                    fName='%s_%sGHzx%sGHz_mean_auto_binned'%(type,f1,f2)
                    lb, cbListAuto, dNameAuto,stdListAuto= increment_dict_list(specDir,fName,'mosaic',[], [] ,[])
                    plot_cl_dict(plotDir,patchName,fName,type,lb,cbListAuto,dNameAuto,fields,theoryFile,lmax,stdList=stdListAuto,Bblfile=None,semilog_T=False,noiseInfo=white_noise_level[f1,f2,'noiseInfo'],beamName=white_noise_level[f1,f2,'beamName'])

                    fName='%s_%sGHzx%sGHz_mean_cross_binned'%(type,f1,f2)
                    lb, cbListCross, dNameCross,stdListCross= increment_dict_list(specDir,fName,'mosaic',[], [],[])
                    plot_cl_dict(plotDir,patchName,fName,type,lb,cbListCross,dNameCross,fields,theoryFile,lmax,stdList=stdListCross,Bblfile=None,semilog_T=False)

                    fName='%s_%sGHzx%sGHz_mean_noise_binned'%(type,f1,f2)
                    lb, NbList, dNameNoise,stdListNoise= increment_dict_list(specDir,fName,'mosaic',[], [],[])
                    plot_cl_dict(plotDir,patchName,fName,type,lb,NbList,dNameNoise,fields,theoryFile,lmax,stdList=stdListNoise,Bblfile=None,semilog_T=True,semilog_pol=True,noiseInfo=white_noise_level[f1,f2,'noiseInfo'],beamName=white_noise_level[f1,f2,'beamName'])

                    count2+=1
                count1+=1
