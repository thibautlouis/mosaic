#!/usr/bin/env python

import healpy as hp
import numpy as np
import pylab as plt
import os

def read_binning_file(file,lmax):
    binLo,binHi,binC = plt.loadtxt(file,unpack=True)
    id = np.where(binHi <lmax)
    binLo,binHi,binC=binLo[id],binHi[id],binC[id]
    if binLo[0]<2:
        binLo[0]=2
    binHi=binHi.astype(np.int)
    binLo=binLo.astype(np.int)
    binSize=binHi-binLo+1
    return (binLo,binHi,binSize)

def get_window_beam_array(p,f1,f2,lmax):
    W={}
    if f1==f2:
        ell,fell_T=np.loadtxt(p['beam_%s_T'%f1],unpack=True)
        ell,fell_Pol=np.loadtxt(p['beam_%s_Pol'%f1],unpack=True)
        ell,fell_T,fell_Pol= ell[:lmax],fell_T[:lmax],fell_Pol[:lmax]
        W['TT']=fell_T**2
        W['TP']=fell_T*fell_Pol
        W['PP']=fell_Pol**2
    else:
        ell,fell_T_1=np.loadtxt(p['beam_%s_T'%f1],unpack=True)
        ell,fell_Pol_1=np.loadtxt(p['beam_%s_Pol'%f1],unpack=True)
        ell,fell_T_2=np.loadtxt(p['beam_%s_T'%f2],unpack=True)
        ell,fell_Pol_2=np.loadtxt(p['beam_%s_Pol'%f2],unpack=True)
        ell,fell_T_1,fell_Pol_1,fell_T_2,fell_Pol_2= ell[:lmax],fell_T_1[:lmax],fell_Pol_1[:lmax],fell_T_2[:lmax],fell_Pol_2[:lmax]
        W['TT']=fell_T_1*fell_T_2
        W['TP']=fell_T_1*fell_Pol_2
        W['PT']=fell_T_2*fell_Pol_1
        W['PP']=fell_Pol_1*fell_Pol_2
    return( ell,W )

def map2alm(map,pixel,lmax=None):
    if pixel=='healpix':
        return hp.sphtfunc.map2alm(map)
    if pixel=='car':
        from enlib import curvedsky
        alm = curvedsky.map2alm(map,lmax= lmax)
        alm = alm.astype(np.complex128)
        return alm

def read_clth_file(clfile,lmax=None):
    fields=['TT','EE','BB','TE']
    cl={}
    lth, cl['TT'], cl['EE'], cl['BB'], cl['TE'] = np.loadtxt(clfile,unpack=True)
    fth=lth*(lth+1)/(2*np.pi)
    for l1 in fields:
        cl[l1]/=fth
        cl[l1][0]=0
    if lmax is None:
        return(lth,cl['TT'], cl['EE'], cl['BB'], cl['TE'])
    else:
        return(lth[:lmax],cl['TT'][:lmax], cl['EE'][:lmax], cl['BB'][:lmax], cl['TE'][:lmax])

def get_cl_dict(alm1,alm2):
    # For the auto it's a bit stupid but here we go
    cl_dict={}
    fields=['TT','EE','BB','TE','EB','TB']
    cls=hp.sphtfunc.alm2cl(alm1,alm2)
    l=np.arange(len(cls[0]))
    count=0
    for l1 in fields:
        cl_dict[l1]=cls[count]
        count+=1
    fields=['TT','EE','BB','ET','BE','BT']
    cls=hp.sphtfunc.alm2cl(alm2,alm1)
    count=0
    for l1 in fields:
        cl_dict[l1]=cls[count]
        count+=1
    return l,cl_dict

def write_cl_dict(dirName,patchName,fName,l,cl_dict,hdf5=False,file=None):
    if hdf5==False:
        fields=['TT','EE','BB','TE','EB','TB','ET','BE','BT']
        cl=[cl_dict[l1] for l1 in fields]
        
        cl[0:0]= [l]
        str='l'
        for l1 in fields:
            str+=' Cl_%s'%l1
        cl=np.array(cl)
        np.savetxt(dirName+'/%s_'%patchName+fName+'.dat', np.transpose(cl),header=str)
    else:
        spec=file.create_group(patchName+fName)
        array=get_cl_array(l,cl_dict)
        spec.create_dataset(name='data',data=array,dtype='float')


def read_cl_dict(dirName,patchName,fName,hdf5=False,file=None):
    
    if hdf5==False:
        data=np.loadtxt(dirName+'/%s_'%patchName+fName+'.dat')
    else:
        spec=file[patchName+'/%s'%fName]
        data=spec['data']

    l=data[:,0]
    cl_dict={}
    fields=['TT','EE','BB','TE','EB','TB','ET','BE','BT']
    count=1
    for l1 in fields:
        cl_dict[l1]=data[:,count]
        count+=1
    return(l,cl_dict)

def bin_cl_dict(l,cl_dict,binningfile,lmax,type):
    
    binLo,binHi,binCent = plt.loadtxt(binningfile,unpack=True)
    id = np.where(binHi<lmax)
    binHi = binHi[id]
    binLo = binLo[id]
    binCent = binCent[id]
    nBins=len(binHi)
    binnedPower={}
    fields=['TT','EE','BB','TE','EB','TB','ET','BE','BT']
    for l1 in fields:
        binnedPower[l1] = binCent.copy()*0.
    
    if type=='Dl':
        f=(l*(l+1)/(2*np.pi))
    if type=='Cl':
        f=l*0+1
    for l1 in fields:
        for ibin in xrange(nBins):
            loc = np.where((l >= binLo[ibin]) & (l <= binHi[ibin]))
            binnedPower[l1][ibin] = (cl_dict[l1][loc]*f[loc]).mean()
    return(binCent,binnedPower)

def mcm_inv_mult(mcm_array,vec,direct_invert=True):
    nb_all=len(vec)
    new_vec=np.zeros(nb_all)
    nb=nb_all/6
    if direct_invert==True:
        vec[:nb]=np.dot(np.linalg.inv(mcm_array[:nb,:nb]),vec[:nb])
        vec[nb:3*nb]=np.dot(np.linalg.inv(mcm_array[nb:3*nb,nb:3*nb]),vec[nb:3*nb])
        vec[3*nb:6*nb]=np.dot(np.linalg.inv(mcm_array[3*nb:6*nb,3*nb:6*nb]),vec[3*nb:6*nb])
    else:
        vec[:nb]=np.linalg.solve(mcm_array[:nb,:nb], vec[:nb])
        vec[nb:3*nb]=np.linalg.solve(mcm_array[nb:3*nb,nb:3*nb], vec[nb:3*nb])
        vec[3*nb:6*nb]=np.linalg.solve(mcm_array[3*nb:6*nb,3*nb:6*nb], vec[3*nb:6*nb])
    return vec


def apply_mcm(cb_dict,mcmlist,direct_invert=False):
    new_cb_dict={}
    nBins=len(cb_dict['TT'])
    fields=[['TT','TE','TB','EE','EB','BB'],['TT','ET','BT','EE','BE','BB']]
    for i in range(2):
        vec=[]
        for l1 in fields[i]:
            vec=np.append(vec,cb_dict[l1])
        
        vec=mcm_inv_mult(mcmlist[i],vec,direct_invert=False)
        count=0
        for l1 in fields[i]:
            new_cb_dict[l1]=vec[count*nBins:(count+1)*nBins]
            count+=1
    return(new_cb_dict)

def get_cl_array(lb,cb_dict):
    array=[]
    array+=[lb]
    fields=['TT','EE','BB','TE','EB','TB','ET','BE','BT']
    for l1 in fields:
        array+=[cb_dict[l1]]
    return(array)

def get_mean_dict(list_dict,get_std=False):
    if len(list_dict) > 1:
        fields=['TT','EE','BB','TE','EB','TB','ET','BE','BT']
        cl_mean={}
        if get_std==True:
            std={}
        for l1 in fields:
            cl_mean[l1]=[]
            for dict in list_dict:
                cl_mean[l1]+=[dict[l1]]
        for l1 in fields:
            if get_std==True:
                std[l1]=np.std(cl_mean[l1],axis=0)
            cl_mean[l1]=np.mean(cl_mean[l1],axis=0)


        if get_std==True:
            return cl_mean,std
        else:
            return cl_mean
    else:
        #        print 'cant take mean of a single dictionnary, return dictionnary'
        return list_dict[0]

def get_noise_dict(auto,cross):
    fields=['TT','EE','BB','TE','EB','TB','ET','BE','BT']
    noise={}
    for l1 in fields:
        noise[l1]=auto[l1]-cross[l1]
    return noise
