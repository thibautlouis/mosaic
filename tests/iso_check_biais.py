#!/usr/bin/env python

import numpy as np
import healpy as hp
import pylab as plt
import iso_ps_utils
import iso_dict
import sys

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
clfile=p['clfile_truc']

lth,clth=iso_ps_utils.read_clth_file(clfile,type=type,lmax=lmax)

auxDir = 'auxMaps_%s/'%pixel
mcmDir= 'mcm_%s/'%pixel
mcDir=  'montecarlo_%s'%pixel

if naMaster==True:
    mcDir +='_namaster'
if removeMean==True:
    mcDir +='_mean_removed'

winList=auxDir+'window_list.txt'
nPatch,freq=iso_window_utils.get_frequency_list(winList)

if hdf5==True:
    file_mc = h5py.File('%s.hdf5'%(mcDir), 'w')
else:
    file_mc= None

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
                                
                        fName='%s_%sGHzx%sGHz_split%dxsplit%d_binned'%(type,f1,f2,s1,s2)
                        meanName='mc_averaged_'+ fName
                        lb,cb_dict=iso_ps_utils.read_cl_dict(mcDir,patchName,meanName,hdf5,file_mc)
                        stdName='std_'+ fName
                        lb,std_dict=iso_ps_utils.read_cl_dict(mcDir,patchName,stdName,hdf5,file_mc)


