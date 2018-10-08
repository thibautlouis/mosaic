#!/usr/bin/env python

import numpy as np
import healpy as hp
import pylab as plt
import iso_ps_utils
import iso_dict
import sys
import iso_window_utils
import iso_map_utils
import iso_spectra_plot_utils
import h5py


def plot_mc_biais(iStart,iStop,plotDir,name,lb,cb_dict,std_dict,clth,cbth,fields):
    
    count=1
    plt.figure(figsize=(28,12))
    for l1 in fields:
        plt.subplot(3,3,count)
        plt.errorbar(lth,clth[l1],label='theory',color='grey')
        plt.errorbar(lb,cbth[l1],label='binned theory',color='black')
        plt.errorbar(lb,cb_dict[l1],std_dict[l1]/np.sqrt(iStop-iStart),fmt='o',label='mean and std')
        plt.ylabel(r'$D^{%s}_{\ell}$'%l1,fontsize=18)
        if count>6:
            plt.xlabel(r'$\ell$',fontsize=22)

        count+=1
    plt.suptitle('%s'%(name),fontsize=22)
    plt.savefig('%s/%s.png'%(plotDir,name),bbox_inches='tight')
    plt.clf()
    plt.close()
    
    count=1
    plt.figure(figsize=(28,12))

    for l1 in fields:
        plt.subplot(3,3,count)
        plt.errorbar(lb,cb_dict[l1]-cbth[l1],std_dict[l1]/np.sqrt(iStop-iStart),fmt='o')
        plt.errorbar(lb,cbth[l1]*0,color='black')
        plt.ylabel(r'$D^{%s, \rm{mean}}_{\ell}- D^{%s, \rm{th}}_{\ell}$'%(l1,l1),fontsize=18)
        if count>6:
            plt.xlabel(r'$\ell$',fontsize=22)
    
        count+=1
    plt.suptitle('%s'%(name),fontsize=22)

    plt.suptitle('diff %s'%(name),fontsize=22)
    plt.savefig('%s/diff_%s.png'%(plotDir,name),bbox_inches='tight')
    plt.clf()
    plt.close()

p = iso_dict.flipperDict()
p.read_from_file(sys.argv[1])

iMin,iMax=p['iMin'],p['iMax']

nSplits= p['nSplit']
pixel=p['pixelisation']
naMaster=p['useNaMaster']
hdf5=p['hdf5']
type=p['type']
removeMean=p['removeMean']
type=p['type']
lmax=p['lmax']
clfile=p['clfile_trunc']
freqTags=p['freqTags']
iStart=p['iStart']
iStop=p['iStop']


auxDir = 'auxMaps_%s/'%pixel
mcmDir= 'mcm_%s/'%pixel
mcDir=  'montecarlo_%s'%pixel
plotDir= 'plots_mc_%s'%pixel

if naMaster==True:
    mcDir +='_namaster'
    plotDir+='_namaster'

if removeMean==True:
    mcDir +='_mean_removed'
    plotDir +='_mean_removed'

iso_map_utils.create_directory(plotDir)


winList=auxDir+'window_list.txt'
nPatch,freq=iso_window_utils.get_frequency_list(winList)

if hdf5==True:
    file_mc = h5py.File('%s.hdf5'%(mcDir), 'r')
else:
    file_mc= None

fields=['TT','EE','BB','TE','EB','TB','ET','BE','BT']

white_noise_level={}
for f1 in freqTags:
    for f2 in freqTags:
        if f1==f2:
            white_noise_level[f1,f2,'noiseInfo']=nSplits,p['rms_%s_T'%f1],p['rms_%s_pol'%f2]
        else :
            white_noise_level[f1,f2]=0,0
        
        white_noise_level[f1,f2,'beamName']=p['beam_%s_T'%f1],p['beam_%s_pol'%f2]


for i in range(nPatch):
    patchName='patch_%03d'%i
    if len(freq[i]) !=0:
        count1=0
        for f1 in freq[i]:
            count2=0
            for f2 in freq[i]:
                if count2>count1:
                    continue
                
                Bbl_T=np.load('mcm_%s/Bbl_T_%03d_%sGHzx%sGHz.npy'%(pixel,i,f1,f2))
                Bbl_Tp=np.load('mcm_%s/Bbl_Tp_%03d_%sGHzx%sGHz.npy'%(pixel,i,f1,f2))
                Bbl_pT=np.load('mcm_%s/Bbl_pT_%03d_%sGHzx%sGHz.npy'%(pixel,i,f1,f2))
                Bbl_pol=np.load('mcm_%s/Bbl_pol_%03d_%sGHzx%sGHz.npy'%(pixel,i,f1,f2))

                count_s1=0
                for s1 in range(nSplits):
                    count_s2=0
                    for s2 in range(nSplits):
                        if ((f1==f2) & (count_s2>count_s1)):
                            break
                        
                        lth,clth=iso_ps_utils.read_clth_file(clfile,type=type,lmax=lmax)
                        cbth=iso_spectra_plot_utils.bin_theory(Bbl_T,Bbl_Tp,Bbl_pT,Bbl_pol,clth,lmax)

                        
                        if (s1==s2) & (f1==f2):
                            nlth,label_th=iso_spectra_plot_utils.get_nlth(white_noise_level[f1,f2,'noiseInfo'],white_noise_level[f1,f2,'beamName'],fields,lmax,type)
                            for l1 in fields:
                                
                                clth[l1]+=nlth[l1]
                            cbth=iso_spectra_plot_utils.bin_theory(Bbl_T,Bbl_Tp,Bbl_pT,Bbl_pol,clth,lmax)
                    
                        fName='%s_%sGHzx%sGHz_split%dxsplit%d_binned'%(type,f1,f2,s1,s2)
                        meanName='mc_averaged_'+ fName
                        lb,cb_dict=iso_ps_utils.read_cl_dict(mcDir,patchName,meanName,hdf5,file_mc)
                        stdName='std_'+ fName
                        lb,std_dict=iso_ps_utils.read_cl_dict(mcDir,patchName,stdName,hdf5,file_mc)
                        name= '%s_%s'%(patchName,meanName)
                        plot_mc_biais(iStart,iStop,plotDir,name,lb,cb_dict,std_dict,clth,cbth,fields)

                        count_s2+= 1
                    count_s1+=1
                count2+=1
            count1+=1
