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
import matplotlib.cm as cm



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

def get_nlth(noiseInfo,beamName,fields,lmax,type):
             
    l,fell_T=np.loadtxt(beamName[0],unpack=True)
    l,fell_Pol=np.loadtxt(beamName[1],unpack=True)
    
    rmsT = noiseInfo[1]
    rmsP = noiseInfo[2]
    nl_th={}
    label_th={}
    for l1 in fields:
        nl_th[l1]=np.zeros(lmax)
        label_th[l1]=''
    
    lth=np.arange(2,lmax+2)

    nl_th['TT']=np.ones(lmax)*(rmsT*np.pi/(60*180))**2/fell_T[:lmax]**2*noiseInfo[0]
    label_th['TT']='rms T %0.2f uk.arcmin'%rmsT
    nl_th['EE']=np.ones(lmax)*(rmsP*np.pi/(60*180))**2/fell_Pol[:lmax]**2*noiseInfo[0]
    nl_th['BB']=np.ones(lmax)*(rmsP*np.pi/(60*180))**2/fell_Pol[:lmax]**2*noiseInfo[0]
    label_th['EE']='rms pol %0.2f uk.arcmin'%rmsP
    label_th['BB']='rms pol %0.2f uk.arcmin'%rmsP
    
    if type=='Dl':
        for l1 in fields:
            nl_th[l1]*=lth*(lth+1)/(2*np.pi)

    return(nl_th,label_th)

def bin_theory(Bbl_T,Bbl_Tp,Bbl_pT,Bbl_pol,clth,lmax):
    cbth={}
    cbth['TT']=np.dot(Bbl_T,clth['TT'])
    cbth['TE']=np.dot(Bbl_Tp,clth['TE'])
    cbth['ET']=np.dot(Bbl_pT,clth['TE'])
    cbth['TB']=cbth['TT'][:lmax]*0
    cbth['BT']=cbth['TT'][:lmax]*0
    
    vec=np.zeros((4*lmax))
    vec[:lmax]=clth['EE'][:lmax]
    vec[lmax:2*lmax]=clth['EE'][:lmax]*0
    vec[2*lmax:3*lmax]=clth['EE'][:lmax]*0
    vec[3*lmax:4*lmax]=clth['BB'][:lmax]
    
    Nbin= Bbl_pol.shape[0]/4
    vecb=np.dot(Bbl_pol,vec)
    cbth['EE']=vecb[:Nbin]
    cbth['EB']=vecb[Nbin:2*Nbin]
    cbth['BE']=vecb[2*Nbin:3*Nbin]
    cbth['BB']=vecb[3*Nbin:4*Nbin]
    return(cbth)


def plot_cl_dict(plotDir,specName,type,lb,dictList,dictNames,fields,theoryFile,lmax,stdList=None,Bblfile=None,\
                 semilog_T=False,semilog_pol=False,noiseInfo=None,beamName=None,cutForPlot=None,all=False,yrange_all=None):


    if noiseInfo is not None and beamName is not None:
        nl_th,label_th=get_nlth(noiseInfo,beamName,fields,lmax,type)

    if theoryFile is not None:
        lth,clth_dict=dict_from_theory_file(fields,theoryFile,lmax)

    if Bblfile is not None:
        print 'Not implemented yet'
        #Bbl=np.load('%s'%Bblfile)
        #cbth_dict=bin_dict(Bbl,clth_dict)

    fb={}

    for l1 in fields:
        if type=='Dl':
            fb[l1]=lb[l1]*0+1
        if type=='Cl':
            fb[l1]=lb[l1]**2/(2*np.pi)

    count=1
    plt.figure(figsize=(24,12))
    for l1 in fields:
        
        idth= np.where(lth< cutForPlot[l1])
        
        if all==False:
            plt.subplot(3,3,count)
        else:
            plt.clf()
            plt.close()
            plt.figure(figsize=(24,12))

                
        if l1=='TT' and semilog_T:
            plt.semilogy()
        if l1=='EE' and semilog_pol:
            plt.semilogy()
        if l1=='BB' and semilog_pol:
            plt.semilogy()
        plt.errorbar(lth[idth],clth_dict[l1][idth],label='CMB',color='black')
        if noiseInfo is not None:
            plt.errorbar(lth[idth],clth_dict[l1][idth]+nl_th[l1][idth],label='CMB+noise',color='blue')
            plt.errorbar(lth[idth],nl_th[l1][idth],label='noise',color='grey')
        if Bblfile:
            plt.errorbar(lb[l1],cbth_dict[l1])
        
        if stdList is not None:
            colorArray = cm.rainbow(np.linspace(0, 1, len(dictList)))
            lshiftArray=np.linspace(-15,15,len(dictList))
            my_count=0
            for dict,std,name in zip(dictList,stdList,dictNames):
                if all==False:
                    plt.errorbar(lb[l1],dict[l1]*fb[l1],std[l1]*fb[l1],fmt='o',label='%s'%name,color='red')
                else:
                    plt.errorbar(lb[l1]+lshiftArray[my_count],dict[l1]*fb[l1],std[l1]*fb[l1],fmt='o',label='%s'%name,color=colorArray[my_count])
                    my_count+=1

        else:
            colorArray = cm.rainbow(np.linspace(0, 1, len(dictList)))
            lshiftArray=np.linspace(-15,15,len(dictList))
            my_count=0
            for dict,name in zip(dictList,dictNames):
                if name=='namaster':
                    plt.errorbar(lb[l1],dict[l1]*fb[l1],fmt='-',label='%s'%name)
                else:
                    if all==False:
                        plt.errorbar(lb[l1],dict[l1]*fb[l1],fmt='o',label='%s'%name,color='red')
                
                    else:
                        plt.errorbar(lb[l1]+lshiftArray[my_count],dict[l1]*fb[l1],fmt='o',label='%s'%name,color=colorArray[my_count])
                        my_count+=1
                            
        plt.ylabel(r'$D^{%s}_{\ell}$'%l1,fontsize=18)

        if all==True:
            if yrange_all is not None:
                plt.ylim(yrange_all[l1][0], yrange_all[l1][1])
            plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            plt.xlabel(r'$\ell$',fontsize=22)
            plt.suptitle('%s_%s_all'%(specName[:9],l1),fontsize=22)
            plt.savefig('%s/%s_%s.png'%(plotDir,specName[:9],l1),bbox_inches='tight')
            plt.clf()
            plt.close()
        else:
            if count==3:
                plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            if count>6:
                plt.xlabel(r'$\ell$',fontsize=22)

            count+=1

    if all==False:

        plt.suptitle('%s'%(specName),fontsize=22)
        plt.savefig('%s/%s.png'%(plotDir,specName),bbox_inches='tight')
        plt.clf()
        plt.close()


def plot_all_spectra(mcmDir,specDir,plotDir,winList,nSplits,hdf5,type,theoryFile,lmax,\
                     compare_mosaic_namaster=False,white_noise_level=None,mcDir=None,cutForPlot=None,yrange_all=None):


    def increment_dict_list(specDir,fName,name, cbList, dictNamesList ,stdList,cutForPlot):

        lb,cb_dict=iso_ps_utils.read_cl_dict(specDir,patchName,fName,hdf5,file)
        lb_field={}
        
        if cutForPlot is not None:
            for l1 in fields:
                id=np.where(lb<cutForPlot[l1])
                cb_dict[l1]=cb_dict[l1][id]
                lb_field[l1]=lb[id]
        else:
            for l1 in fields:
                lb_field[l1]=lb
    
        cbList+=[cb_dict]
        dictNamesList+=[name]
        if mcDir is not None:
            lb,std_dict=iso_ps_utils.read_cl_dict(mcDir,patchName,'std_'+fName,hdf5,file_mc)
            if cutForPlot is not None:
                for l1 in fields:
                    id=np.where(lb<cutForPlot[l1])
                    std_dict[l1]=std_dict[l1][id]


            stdList+=[std_dict]
        else:
            stdList=None
        
        
        return (lb_field, cbList, dictNamesList,stdList)


    fields=['TT','EE','BB','TE','EB','TB','ET','BE','BT']

    if hdf5==True:
        file = h5py.File('%s.hdf5'%(specDir), 'r')
        if compare_mosaic_namaster==True:
            file_nm = h5py.File('%s'%(specDir)+'_namaster.hdf5', 'r')
        if mcDir is not None:
            file_mc=h5py.File('%s.hdf5'%(mcDir), 'r')
    else:
        file=None
        file_nm=None

    nPatch,freq=iso_window_utils.get_frequency_list(winList)
    for i in range(nPatch):
        
        patchName='patch_%03d'%i
        
        cbList_all=[]
        dictNamesList_all=[]
        if mcDir is not None:
            stdList_all=[]
        else:
            stdList_all=None

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
                            lb, cbList, dictNamesList,stdList= increment_dict_list(specDir,fName,'s%d %sx s%d %s'%(s1,f1,s2,f2),cbList, dictNamesList ,stdList,cutForPlot)
                            if not (f1==f2) & (count_s2==count_s1):
                                lb, cbList_all, dictNamesList_all,stdList_all= increment_dict_list(specDir,fName,'s%d %sx s%d %s'%(s1,f1,s2,f2),cbList_all, dictNamesList_all ,stdList_all,cutForPlot)

                            if compare_mosaic_namaster:
                                id= specDir[-4:]
                                naMasterSpecDir=specDir[:-4]+'_namaster'+id
                                lb, cbList, dictNamesList,stdList=increment_dict_list(naMasterSpecDir,fName,'namaster',cbList, dictNamesList ,stdList,cutForPlot)

                            specName= '%s_%s'%(patchName,fName)
                            if (f1==f2) & (s1==s2):
                                plot_cl_dict(plotDir,specName,type,lb,cbList,dictNamesList,fields,theoryFile,lmax,stdList=stdList,Bblfile=None,\
                                             semilog_T=True,noiseInfo=white_noise_level[f1,f2,'noiseInfo'],beamName=white_noise_level[f1,f2,'beamName'],cutForPlot=cutForPlot)

                            else:
                                plot_cl_dict(plotDir,specName,type,lb,cbList,dictNamesList,fields,theoryFile,lmax,stdList=stdList,Bblfile=None,\
                                             semilog_T=True,cutForPlot=cutForPlot)
          
                            count_s2+=1
                        count_s1+=1
                        
                    fName='%s_%sGHzx%sGHz_mean_auto_binned'%(type,f1,f2)
                    specName= '%s_%s'%(patchName,fName)

                    lb, cbListAuto, dNameAuto,stdListAuto= increment_dict_list(specDir,fName,'mosaic',[], [] ,[],cutForPlot)
                    plot_cl_dict(plotDir,specName,type,lb,cbListAuto,dNameAuto,fields,theoryFile,lmax,stdList=stdListAuto,Bblfile=None,\
                                 semilog_T=False,noiseInfo=white_noise_level[f1,f2,'noiseInfo'],beamName=white_noise_level[f1,f2,'beamName'],cutForPlot=cutForPlot)

                    fName='%s_%sGHzx%sGHz_mean_cross_binned'%(type,f1,f2)
                    specName= '%s_%s'%(patchName,fName)

                    lb, cbListCross, dNameCross,stdListCross= increment_dict_list(specDir,fName,'mosaic',[], [],[],cutForPlot)
                    plot_cl_dict(plotDir,specName,type,lb,cbListCross,dNameCross,fields,theoryFile,lmax,stdList=stdListCross,Bblfile=None,\
                                 semilog_T=False,cutForPlot=cutForPlot)

                    fName='%s_%sGHzx%sGHz_mean_noise_binned'%(type,f1,f2)
                    specName= '%s_%s'%(patchName,fName)
                    lb, NbList, dNameNoise,stdListNoise= increment_dict_list(specDir,fName,'mosaic',[], [],[],cutForPlot)
                    
                    if f1==f2:
                        plot_cl_dict(plotDir,specName,type,lb,NbList,dNameNoise,fields,theoryFile,lmax,stdList=stdListNoise,Bblfile=None,\
                                     semilog_T=True,semilog_pol=True,noiseInfo=white_noise_level[f1,f2,'noiseInfo'],beamName=white_noise_level[f1,f2,'beamName'],cutForPlot=cutForPlot)
                    else:
                        plot_cl_dict(plotDir,specName,type,lb,NbList,dNameNoise,fields,theoryFile,lmax,stdList=stdListNoise,Bblfile=None,\
                                     semilog_T=False,semilog_pol=False,noiseInfo=white_noise_level[f1,f2,'noiseInfo'],beamName=white_noise_level[f1,f2,'beamName'],cutForPlot=cutForPlot)

                    count2+=1
                count1+=1

            plot_cl_dict(plotDir,specName,type,lb,cbList_all,dictNamesList_all,fields,theoryFile,lmax,stdList=stdList_all,Bblfile=None,semilog_T=False,cutForPlot=cutForPlot,all=True,yrange_all=yrange_all)

