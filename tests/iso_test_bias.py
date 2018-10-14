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
import iso_map_plot_utils
import iso_html_utils
import h5py
import os
from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw

def plot_mc_biais(iStart,iStop,plotDir,name,lb,cb_dict,std_dict,lth,clth,cbth,fields,semilogy_T=False,semilogy_pol=False,cutForPlot=None):
    
    lth_field={}
    lb_field={}

    if cutForPlot is not None:
        for l1 in fields:
            id=np.where(lb<cutForPlot[l1])
            lb_field[l1]=lb[id]
            lth_field[l1]=lth[:cutForPlot[l1]]
            clth[l1]=clth[l1][:cutForPlot[l1]]
            cbth[l1]=cbth[l1][id]
            cb_dict[l1]=cb_dict[l1][id]
            std_dict[l1]=std_dict[l1][id]
    else:
        for l1 in fields:
            lth_field[l1]=lth
            lb_field[l1]=lb


    count=1
    plt.figure(figsize=(28,12))
    for l1 in fields:
        plt.subplot(3,3,count)
        if semilogy_T==True and l1=='TT':
            plt.semilogy()
        if semilogy_pol==True and l1=='EE':
            plt.semilogy()
        if semilogy_pol==True and l1=='BB':
            plt.semilogy()

        
        plt.errorbar(lth_field[l1],clth[l1],label='theory',color='grey')
        plt.errorbar(lb_field[l1],cbth[l1],label='binned theory',color='black')
        plt.errorbar(lb_field[l1],cb_dict[l1],std_dict[l1],fmt='o',label='mean and std')
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
        plt.errorbar(lb_field[l1],cb_dict[l1]-cbth[l1],std_dict[l1]/np.sqrt(iStop-iStart),fmt='o')
        plt.errorbar(lb_field[l1],cbth[l1]*0,color='black')
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
mask=p['mask']
survey_mask_coordinates=p['survey_mask_coordinate']
tessel_healpix=p['tessel_healpix']
noise=p['noise']

mapDir= 'maps_%s'%pixel
auxDir = 'auxMaps_%s/'%pixel
mcmDir= 'mcm_%s/'%pixel
mcDir=  'montecarlo_%s'%pixel
plotDir= 'plots_mc_%s'%pixel
htmlDir= 'webpage_mc_%s'%pixel

if naMaster==True:
    mcDir +='_namaster'
    plotDir+='_namaster'
    htmlDir+='_namaster'
if removeMean==True:
    mcDir +='_mean_removed'
    plotDir +='_mean_removed'
    htmlDir +='_mean_removed'
mapDir+='_%03d/'%(0)

iso_map_utils.create_directory(plotDir)
iso_map_utils.create_directory(htmlDir)


winList=auxDir+'window_list.txt'
#iso_map_plot_utils.plot_survey_map(auxDir,mapDir,plotDir,pixel,winList,mask,freqTags,survey_mask_coordinates=survey_mask_coordinates,tessel_healpix=tessel_healpix,color_range=None)
#iso_map_plot_utils.plot_all_windows(auxDir,plotDir,pixel,winList,tessel_healpix=tessel_healpix)

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
            white_noise_level[f1,f2,'noiseInfo']=nSplits,noise['rms_%s_T'%f1],noise['rms_%s_pol'%f2]
        else :
            white_noise_level[f1,f2]=0,0
        
        white_noise_level[f1,f2,'beamName']=p['beam_%s_T'%f1],p['beam_%s_pol'%f2]

cutForPlot={}
cutForPlot['TT']=2000
cutForPlot['EE']=1400
cutForPlot['TE']=1500
cutForPlot['ET']=1500
cutForPlot['BB']=1400
cutForPlot['EB']=1400
cutForPlot['BE']=1400
cutForPlot['BT']=1400
cutForPlot['TB']=1400


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
                            #print white_noise_level[f1,f2,'noiseInfo'],white_noise_level[f1,f2,'beamName'],fields,lmax,type
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
                        if s1==s2:
                            plot_mc_biais(iStart,iStop,plotDir,name,lb,cb_dict,std_dict,lth,clth,cbth,fields,semilogy_T=True,semilogy_pol=True,cutForPlot=cutForPlot)
                        else:
                            plot_mc_biais(iStart,iStop,plotDir,name,lb,cb_dict,std_dict,lth,clth,cbth,fields,semilogy_T=False,semilogy_pol=False,cutForPlot=cutForPlot)

                        count_s2+= 1
                    count_s1+=1
                count2+=1
            count1+=1


path=(os.path.dirname(os.path.realpath(__file__)))
os.system('cp %s/../data/multistep2.js ./%s/multistep2.js'%(path,htmlDir))
    
fileName='%s/isotropy.html'%htmlDir
g = open(fileName,mode="w")
g.write('<html>\n')
g.write('<head>\n')
g.write('<title> Mosaic </title>\n')
g.write('<script src="multistep2.js"></script>\n')
g.write('<script> add_step("sub",  ["j","k"]) </script> \n')
g.write('<script> add_step("all",  ["c","v"]) </script> \n')
g.write('<script> add_step("diff",  ["a","z"]) </script> \n')
g.write('</head> \n')
g.write('<body> \n')
g.write('<div class=sub> \n')
    
for i in range(nPatch):
    patchName='patch_%03d'%i
    g.write('<div class=all>\n')
    
    survey= Image.open('%s/im_%03d.png'%(plotDir,i))
    count1=0
    for f1 in freqTags:
        count2=0
        for f2 in freqTags:
            if count2>count1:
                continue
            count_s1=0
            for s1 in range(nSplits):
                count_s2=0
                for s2 in range(nSplits):
                    if ((f1==f2) & (count_s2>count_s1)):
                        break
                    
                    g.write('<div class=diff> \n')

                    fName='patch_%03d_mc_averaged_%s_%sGHzx%sGHz_split%dxsplit%d_binned.png'%(i,type,f1,f2,s1,s2)
                    iso_html_utils.spectra_image(htmlDir,plotDir,fName,i,null=True)
                    g.write('<img src="'+fName+'" width="100%" /> \n')
                    fName='diff_patch_%03d_mc_averaged_%s_%sGHzx%sGHz_split%dxsplit%d_binned.png'%(i,type,f1,f2,s1,s2)
                    iso_html_utils.spectra_image(htmlDir,plotDir,fName,i,null=True)
                    g.write('<img src="'+fName+'" width="100%" /> \n')
                    
                    g.write('</div>\n')


                    count_s2+=1
                count_s1+=1
            count2+=1
        count1+=1
        
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
                        
                        g.write('<div class=diff> \n')

                        fName='patch_%03d_mc_averaged_%s_%sGHzx%sGHz_split%dxsplit%d_binned.png'%(i,type,f1,f2,s1,s2)
                        iso_html_utils.spectra_image(htmlDir,plotDir,fName,i,null=False)
                        g.write('<img src="'+fName+'" width="100%" /> \n')
                        fName='diff_patch_%03d_mc_averaged_%s_%sGHzx%sGHz_split%dxsplit%d_binned.png'%(i,type,f1,f2,s1,s2)
                        iso_html_utils.spectra_image(htmlDir,plotDir,fName,i,null=False)
                        g.write('<img src="'+fName+'" width="100%" /> \n')

                        g.write('</div>\n')

                        count_s2+=1
                    count_s1+=1
                count2+=1
            count1+=1
    
    g.write('</div>\n')

g.write('</div> \n')
g.write('</body> \n')
g.write('</html> \n')
g.close()
