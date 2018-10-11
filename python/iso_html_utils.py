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
from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw

im_size_healpix=(2200,1900)

def spectra_image(htmlDir,plotDir,fName,i,null=False):
    
    path=(os.path.dirname(os.path.realpath(__file__)))
    full_im = Image.new("RGB", im_size_healpix,color='white')
    survey= Image.open('%s/im_%03d.png'%(plotDir,i))
    w_survey, h_survey = survey.size
    x0_survey,y0_survey=538,100
    full_im.paste(survey, (x0_survey,y0_survey,x0_survey+w_survey,y0_survey+h_survey))
    draw = ImageDraw.Draw(full_im)
    font = ImageFont.truetype('%s/../data/arial.ttf'%path, 48, encoding="unic")

    
    if null==True:
        text='No %s spectra available for this patch'%fName
        draw.text((0, 800),text,(0,0,255),font)
        
        text='press j/k to go from one patch to the other, c/v to go from one data product to the other '
        draw.text((170, 0),text,(0,0,0),font)

        full_im.save(os.path.expanduser('%s/%s'%(htmlDir,fName)))
    else:
        text=' press j/k to go from one patch to the other, c/v to go from one data product to the other '
        draw.text((170, 0),text,(0,0,0),font)
        spec= Image.open('%s/%s'%(plotDir,fName))
        w_spec, h_spec = spec.size
        x0_spec,y0_spec=0,h_survey+y0_survey
        full_im.paste(spec, (x0_spec,y0_spec,x0_spec+w_spec,y0_spec+h_spec))
        full_im.save(os.path.expanduser('%s/%s'%(htmlDir,fName)))


def make_html_healpix(htmlDir,plotDir,freqTags,winList,nSplits,type):
    
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
    g.write('</head> \n')
    g.write('<body> \n')
    g.write('<div class=sub> \n')

    nPatch,freq=iso_window_utils.get_frequency_list(winList)
    
    for i in range(nPatch):
        
        patchName='patch_%03d'%i
        
        g.write('<div class=all>\n')

        count1=0
        for f1 in freqTags:
            full_im = Image.new("RGB", im_size_healpix,color='white')
            survey= Image.open('%s/im_%03d.png'%(plotDir,i))
            w_survey, h_survey = survey.size
            x0_survey,y0_survey=538,100
            full_im.paste(survey, (x0_survey,y0_survey,x0_survey+w_survey,y0_survey+h_survey))
            draw = ImageDraw.Draw(full_im)
            font = ImageFont.truetype('%s/../data/arial.ttf'%path, 48, encoding="unic")
            text='press j/k to go from one patch to the other, c/v to go from one data product to the other '
            draw.text((170, 0),text,(0,0,0),font)

            
            text='No %s GHz data available for this patch (due to mask)'%f1
            draw.text((300, 800),text,(0,0,255),font)
            str='image_%03d_%s.png'%(i,f1)
            full_im.save(os.path.expanduser('%s/%s'%(htmlDir,str)))
            g.write('<img src="'+str+'" width="100%" /> \n')
            
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
                        fName='patch_%03d_%s_%sGHzx%sGHz_split%dxsplit%d_binned.png'%(i,type,f1,f2,s1,s2)
                        spectra_image(htmlDir,plotDir,fName,i,null=True)
                        g.write('<img src="'+fName+'" width="100%" /> \n')

                        count_s2+=1
                    count_s1+=1
                        
                fName='patch_%03d_%s_%sGHzx%sGHz_mean_auto_binned.png'%(i,type,f1,f2)
                spectra_image(htmlDir,plotDir,fName,i,null=True)
                g.write('<img src="'+fName+'" width="100%" /> \n')

                fName='patch_%03d_%s_%sGHzx%sGHz_mean_cross_binned.png'%(i,type,f1,f2)
                spectra_image(htmlDir,plotDir,fName,i,null=True)
                g.write('<img src="'+fName+'" width="100%" /> \n')

                fName='patch_%03d_%s_%sGHzx%sGHz_mean_noise_binned.png'%(i,type,f1,f2)
                spectra_image(htmlDir,plotDir,fName,i,null=True)
                g.write('<img src="'+fName+'" width="100%" /> \n')


                count2+=1
            count1+=1

        if len(freq[i]) !=0:
            
            count1=0
            for f1 in freq[i]:
                winName_T= 'window_%s_%03d_T_projected'%(f1,i)
                winName_pol= 'window_%s_%03d_pol_projected'%(f1,i)

                full_im = Image.new("RGB", im_size_healpix,color='white')
                
                win_T= Image.open('%s/%s.png'%(plotDir,winName_T))
                w_win, h_win = win_T.size
                x0_win_T,y0_win_T=0,100
                full_im.paste(win_T, (x0_win_T,y0_win_T,x0_win_T+w_win,y0_win_T+h_win))

                survey= Image.open('%s/im_%03d.png'%(plotDir,i))
                w_survey, h_survey = survey.size
                x0_survey,y0_survey=w_win,100
                full_im.paste(survey, (x0_survey,y0_survey,x0_survey+w_survey,y0_survey+h_survey))
                
                win_Pol= Image.open('%s/%s.png'%(plotDir,winName_pol))
                x0_win_pol,y0_win_pol=w_survey+w_win,100
                full_im.paste(win_Pol,(x0_win_pol,y0_win_pol,x0_win_pol+w_win,y0_win_pol+h_win))

                for s in range(nSplits):
                    y0=y0_win_pol+(s+1)*h_win
                    count=0
                    for l1 in ['T','Q','U']:
                        x0=100+w_win*count
                        im=Image.open('%s/split_%d_%s_%03d_%s_projected.png'%(plotDir,s,f1,i,l1))
                        print im.size
                        print w_win,h_win
                        coord=(x0, y0,x0+w_win, y0+h_win)

                        
                        full_im.paste(im,coord)
                        count+=1
                draw = ImageDraw.Draw(full_im)
                text='press j/k to go from one patch to the other, c/v to go from one data product to the other '
                draw.text((170, 0),text,(0,0,0),font)

                str='image_%03d_%s.png'%(i,f1)
                full_im.save(os.path.expanduser('%s/%s'%(htmlDir,str)))
        
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
                            fName='patch_%03d_%s_%sGHzx%sGHz_split%dxsplit%d_binned.png'%(i,type,f1,f2,s1,s2)
                            spectra_image(htmlDir,plotDir,fName,i)

                            count_s2+=1
                        count_s1+=1
                                                            
                    fName='patch_%03d_%s_%sGHzx%sGHz_mean_auto_binned.png'%(i,type,f1,f2)
                    spectra_image(htmlDir,plotDir,fName,i)
                        
                    fName='patch_%03d_%s_%sGHzx%sGHz_mean_cross_binned.png'%(i,type,f1,f2)
                    spectra_image(htmlDir,plotDir,fName,i)
                        
                    fName='patch_%03d_%s_%sGHzx%sGHz_mean_noise_binned.png'%(i,type,f1,f2)
                    spectra_image(htmlDir,plotDir,fName,i)
                    
                    count2+=1
                count1+=1

        g.write('</div>\n')
        
    g.write('</div> \n')
    g.write('</body> \n')
    g.write('</html> \n')
    g.close()
