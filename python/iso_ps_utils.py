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
    ell,fell_T_1=np.loadtxt(p['beam_%s_T'%f1],unpack=True)
    ell,fell_Pol_1=np.loadtxt(p['beam_%s_pol'%f1],unpack=True)
    ell,fell_T_2=np.loadtxt(p['beam_%s_T'%f2],unpack=True)
    ell,fell_Pol_2=np.loadtxt(p['beam_%s_pol'%f2],unpack=True)
    ell,fell_T_1,fell_Pol_1,fell_T_2,fell_Pol_2= ell[:lmax],fell_T_1[:lmax],fell_Pol_1[:lmax],fell_T_2[:lmax],fell_Pol_2[:lmax]
    W['TT']=fell_T_1*fell_T_2
    W['TP']=fell_T_1*fell_Pol_2
    W['PT']=fell_T_2*fell_Pol_1
    W['PP']=fell_Pol_1*fell_Pol_2
    
    return( ell,W )

def map2alm(map,pixel,niter,lmax=None,theta_range=None):
    
    # This is a generalized version of map2alm"
    # Allow to do spherical transform on healpix, healpix with theta_cut and car
    # Also allow to iterate to get more precise transform (this bit need more testing)
    
    if pixel=='healpix':
        if theta_range is None:
            alm= hp.sphtfunc.map2alm(map,iter=niter)
        else:
            from enlib import curvedsky
            nside=hp.pixelfunc.get_nside(map)
            alm= curvedsky.map2alm_healpix(map,lmax=3*nside-1,theta_min=theta_range[0], theta_max=theta_range[1])
            if iter !=0:
                map_copy=map.copy()
                alm= curvedsky.map2alm_healpix(map,lmax=3*nside-1,theta_min=theta_range[0], theta_max=theta_range[1])
                for k in range(niter):
                    alm += curvedsky.map2alm_healpix(map-curvedsky.alm2map_healpix(alm,map_copy),lmax=3*nside-1,theta_min=thetas[0], theta_max=thetas[1])
        return alm

    if pixel=='car':
        from enlib import curvedsky
        alm = curvedsky.map2alm(map,lmax= lmax)
        if iter !=0:
            map_copy=map.copy()
            for k in range(niter):
                alm += curvedsky.map2alm(map-curvedsky.alm2map(alm,map_copy),lmax=lmax)

        alm = alm.astype(np.complex128)
        return alm

def read_clth_file(clfile,type,lmax=None):
    fields=['TT','EE','BB','TE']
    cl={}
    lth, cl['TT'], cl['EE'], cl['BB'], cl['TE'] = np.loadtxt(clfile,unpack=True)
    fth=lth*(lth+1)/(2*np.pi)
    if type=='Cl':
        for l1 in fields:
            cl[l1]/=fth
            cl[l1][0]=0

    if lmax is not None:
        lth=lth[:lmax]
        for l1 in fields:
            cl[l1]=cl[l1][:lmax]

    cl['ET']=cl['TE']
    cl['TB']=cl['TT']*0
    cl['BT']=cl['TT']*0
    cl['EB']=cl['TT']*0
    cl['BE']=cl['TT']*0

    return(lth,cl)

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
        spec=file.create_group(patchName+'_'+fName)
        array=get_cl_array(l,cl_dict)
        spec.create_dataset(name='data',data=array,dtype='float')


def read_cl_dict(dirName,patchName,fName,hdf5=False,file=None):
    
    if hdf5==False:
        data=np.loadtxt(dirName+'/%s_'%patchName+fName+'.dat')
    else:
        spec=file[patchName+'_'+fName]
        data=np.array(spec['data']).T

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
    nb=nb_all/9
    if direct_invert==True:
        vec[:nb]=np.dot(np.linalg.inv(mcm_array[:nb,:nb]),vec[:nb])
        vec[nb:5*nb]=np.dot(np.linalg.inv(mcm_array[nb:5*nb,nb:5*nb]),vec[nb:5*nb])
        vec[5*nb:9*nb]=np.dot(np.linalg.inv(mcm_array[5*nb:9*nb,5*nb:9*nb]),vec[5*nb:9*nb])
    else:
        vec[:nb]=np.linalg.solve(mcm_array[:nb,:nb], vec[:nb])
        vec[nb:5*nb]=np.linalg.solve(mcm_array[nb:5*nb,nb:5*nb], vec[nb:5*nb])
        vec[5*nb:9*nb]=np.linalg.solve(mcm_array[5*nb:9*nb,5*nb:9*nb], vec[5*nb:9*nb])
    return vec


def apply_mcm(cb_dict,mcm,direct_invert=False):
    new_cb_dict={}
    nBins=len(cb_dict['TT'])
    fields=['TT','TE','TB','ET','BT','EE','EB','BE','BB']
    vec=[]
    for l1 in fields:
        vec=np.append(vec,cb_dict[l1])
        
    vec=mcm_inv_mult(mcm,vec,direct_invert=False)
    count=0
    for l1 in fields:
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
