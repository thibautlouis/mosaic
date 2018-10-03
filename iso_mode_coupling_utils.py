import healpy as hp
import numpy as np
import pylab as plt
import iso_map_utils
import iso_ps_utils
import iso_dict
import sys
import os
import mcm_code_full


def bin_MCM(mcm,mcm_Tp,mcm_pT,mcm_pp,mcm_mm,binLo,binHi,binsize,nBins,scalefac):
    mbb=np.zeros((nBins,nBins))
    mbb_Tp=np.zeros((nBins,nBins))
    mbb_pT=np.zeros((nBins,nBins))
    mbb_pp=np.zeros((nBins,nBins))
    mbb_mm=np.zeros((nBins,nBins))
    mcm_code_full.bin_mcm(mcm.T, binLo,binHi,binsize.T, mbb.T,scalefac)
    mcm_code_full.bin_mcm(mcm_Tp.T, binLo,binHi,binsize.T, mbb_Tp.T,scalefac)
    mcm_code_full.bin_mcm(mcm_pT.T, binLo,binHi,binsize.T, mbb_pT.T,scalefac)
    mcm_code_full.bin_mcm(mcm_pp.T, binLo,binHi,binsize.T, mbb_pp.T,scalefac)
    mcm_code_full.bin_mcm(mcm_mm.T, binLo,binHi,binsize.T, mbb_mm.T,scalefac)
    return(mbb,mbb_Tp,mbb_pT,mbb_pp,mbb_mm)
    
def save_mbb(mcmDir,i,f1,f2,mbb,mbb_Tp,mbb_pT,mbb_pp,mbb_mm,nBins,extra='normal'):
    
    mbb_all=np.zeros((9*nBins,9*nBins))
    mbb_all[:nBins,:nBins]=mbb
    mbb_all[nBins:2*nBins,nBins:2*nBins]=mbb_Tp
    mbb_all[2*nBins:3*nBins,2*nBins:3*nBins]=mbb_Tp
    mbb_all[3*nBins:4*nBins,3*nBins:4*nBins]=mbb_pT
    mbb_all[4*nBins:5*nBins,4*nBins:5*nBins]=mbb_pT
    mbb_all[5*nBins:6*nBins,5*nBins:6*nBins]=mbb_pp
    mbb_all[6*nBins:7*nBins,6*nBins:7*nBins]=mbb_pp
    mbb_all[6*nBins:7*nBins,7*nBins:8*nBins]=-mbb_mm
    mbb_all[7*nBins:8*nBins,6*nBins:7*nBins]=-mbb_mm
    mbb_all[7*nBins:8*nBins,7*nBins:8*nBins]=mbb_pp
    mbb_all[8*nBins:9*nBins,5*nBins:6*nBins]=mbb_mm
    mbb_all[5*nBins:6*nBins,8*nBins:9*nBins]=mbb_mm
    mbb_all[8*nBins:9*nBins,8*nBins:9*nBins]=mbb_pp

    
    np.savetxt('%s/mode_coupling_%03d_%sGHzx%sGHz.dat'%(mcmDir,i,f1,f2),mbb_all)
    return mbb_all



def get_BBl(mcmDir,i,f1,f2, mcm,mcm_Tp,mcm_pT,mcm_pp,mcm_mm,mbb_all,nBins,lmax,binLo,binHi,binSize,extra='normal',scalefac=0):
    
    def inv_mbb_pol_array(mbb_all):
        mbb_pol=mbb_all[5*nBins:9*nBins,5*nBins:9*nBins]
        inv=np.linalg.inv(mbb_pol)
        return inv
    
    def Bbl_pol_array(Bbl_pp,Bbl_mm):
        
        Bbl_pol=np.zeros((4*nBins,4*lmax))
        Bbl_pol[:nBins,:lmax]= Bbl_pp
        Bbl_pol[3*nBins:4*nBins,3*lmax:4*lmax]= Bbl_pp
        Bbl_pol[3*nBins:4*nBins,:lmax]= Bbl_mm
        Bbl_pol[:nBins,3*lmax:4*lmax]= Bbl_mm
        
        Bbl_pol[nBins:2*nBins,lmax:2*lmax]= Bbl_pp
        Bbl_pol[nBins:2*nBins,3*lmax:4*lmax]= -Bbl_mm
        Bbl_pol[2*nBins:3*nBins,lmax:2*lmax]= -Bbl_mm
        Bbl_pol[2*nBins:3*nBins,2*lmax:3*lmax]= Bbl_pp

        return(Bbl_pol)
                
    Bbl=np.zeros((nBins,lmax))
    Bbl_pT=np.zeros((nBins,lmax))
    Bbl_Tp=np.zeros((nBins,lmax))
    Bbl_pp=np.zeros((nBins,lmax))
    Bbl_mm=np.zeros((nBins,lmax))

    mcm_code_full.binning_matrix(mcm.T,binLo,binHi,binSize, Bbl.T,scalefac)
    mcm_code_full.binning_matrix(mcm_Tp.T,binLo,binHi,binSize, Bbl_Tp.T,scalefac)

    mcm_code_full.binning_matrix(mcm_pT.T,binLo,binHi,binSize, Bbl_pT.T,scalefac)

    mcm_code_full.binning_matrix(mcm_pp.T,binLo,binHi,binSize, Bbl_pp.T,scalefac)
    mcm_code_full.binning_matrix(mcm_mm.T,binLo,binHi,binSize, Bbl_mm.T,scalefac)
    
    Bbl_T=np.dot(np.linalg.inv(mbb_all[:nBins,:nBins]),Bbl)
    Bbl_Tp=np.dot(np.linalg.inv(mbb_all[nBins:2*nBins,nBins:2*nBins]),Bbl_Tp)
    Bbl_pT=np.dot(np.linalg.inv(mbb_all[3*nBins:4*nBins,3*nBins:4*nBins]),Bbl_pT)

    Bbl_pol= Bbl_pol_array(Bbl_pp,Bbl_mm)
    inv_pol= inv_mbb_pol_array(mbb_all)
    Bbl_pol=np.dot(inv_pol,Bbl_pol)
    
    np.save('%s/Bbl_T_%03d_%sGHzx%sGHz'%(mcmDir,i,f1,f2),Bbl_T)
    np.save('%s/Bbl_Tp_%03d_%sGHzx%sGHz'%(mcmDir,i,f1,f2),Bbl_Tp)
    np.save('%s/Bbl_pT_%03d_%sGHzx%sGHz'%(mcmDir,i,f1,f2),Bbl_pT)
    np.save('%s/Bbl_pol_%03d_%sGHzx%sGHz'%(mcmDir,i,f1,f2),Bbl_pol)




def get_mode_coupling(auxDir,mcmDir,i,f1,f2,Wl_array,binLo,binHi,binSize,lmax,type,pixel):
    
    nBins=len(binHi)
    
    def init_mcm(lmax):
        mcm=np.zeros((lmax,lmax))
        mcm_Tp=np.zeros((lmax,lmax))
        mcm_pT=np.zeros((lmax,lmax))
        mcm_pp=np.zeros((lmax,lmax))
        mcm_mm=np.zeros((lmax,lmax))
        return(mcm,mcm_Tp,mcm_pT,mcm_pp,mcm_mm)
    
    if type=='Dl':
        scalefac=1
    if type=='Cl':
        scalefac=0


    fields=['TT','TP','PT','PP']

    win_T_1=iso_map_utils.read_map(auxDir,'/window_%s_%03d_T.fits'%(f1,i),pixel)
    win_pol_1=iso_map_utils.read_map(auxDir,'/window_%s_%03d_pol.fits'%(f1,i),pixel)
    win_T_2=iso_map_utils.read_map(auxDir,'/window_%s_%03d_T.fits'%(f2,i),pixel)
    win_pol_2=iso_map_utils.read_map(auxDir,'/window_%s_%03d_pol.fits'%(f2,i),pixel)

    alm_T_1=iso_ps_utils.map2alm(win_T_1,pixel,niter=3,lmax=lmax)
    alm_Pol_1=iso_ps_utils.map2alm(win_pol_1,pixel,niter=3,lmax=lmax)
    alm_T_2=iso_ps_utils.map2alm(win_T_2,pixel,niter=3,lmax=lmax)
    alm_Pol_2=iso_ps_utils.map2alm(win_pol_2,pixel,niter=3,lmax=lmax)
        
    wcl={}
        
    wcl['TT']= hp.alm2cl(alm_T_1,alm_T_2)
    wcl['TP']= hp.alm2cl(alm_T_1,alm_Pol_2)
    wcl['PT']= hp.alm2cl(alm_T_2,alm_Pol_1)
    wcl['PP']= hp.alm2cl(alm_Pol_1,alm_Pol_2)

    l=np.arange(len(wcl['TT']))

    for l1 in fields:
        wcl[l1]*=(2*l+1)
        
    mcm,mcm_pT,mcm_Tp,mcm_pp,mcm_mm=init_mcm(lmax)
    mcm_code_full.calc_mcm(wcl['TT'],wcl['TP'],wcl['PT'],wcl['PP'], Wl_array['TT'],Wl_array['TP'],Wl_array['PT'], Wl_array['PP'],mcm.T,mcm_Tp.T,mcm_pT.T,mcm_pp.T,mcm_mm.T)
    mbb,mbb_Tp,mbb_pT,mbb_pp,mbb_mm=bin_MCM(mcm,mcm_Tp,mcm_pT,mcm_pp,mcm_mm,binLo,binHi,binSize,nBins,scalefac)
    mbb_all=save_mbb(mcmDir,i,f1,f2,mbb,mbb_Tp,mbb_pT,mbb_pp,mbb_mm,nBins)
    get_BBl(mcmDir,i,f1,f2, mcm,mcm_Tp,mcm_pT,mcm_pp,mcm_mm,mbb_all,nBins,lmax,binLo,binHi,binSize,scalefac=scalefac)



