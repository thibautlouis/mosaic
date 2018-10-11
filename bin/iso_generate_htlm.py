#!/usr/bin/env python


import healpy as hp
import numpy as np
import pylab as plt
import iso_dict
import iso_html_utils
import iso_map_utils
import sys
import h5py
import os

p = iso_dict.flipperDict()
p.read_from_file(sys.argv[1])

freqTags=p['freqTags']
nSplits= p['nSplit']
type= p['type']
pixel=p['pixelisation']
removeMean=p['removeMean']
naMaster=p['useNaMaster']
plotDir= 'plots_%s'%pixel
htmlDir= 'webpage_%s'%pixel

auxDir = 'auxMaps_%s/'%pixel

if naMaster==True:
    plotDir+='_namaster'
    htmlDir+='_namaster'
if removeMean==True:
    plotDir +='_mean_removed'
    htmlDir+='_mean_removed'
if len(sys.argv)> 2:
    plotDir += '_%s'%(sys.argv[2])
    htmlDir += '_%s'%(sys.argv[2])


iso_map_utils.create_directory(htmlDir)

winList=auxDir+'window_list.txt'

if pixel=='healpix':
    iso_html_utils.make_html_healpix(htmlDir,plotDir,freqTags,winList,nSplits,type)
if pixel=='car':
    print 'not implemented yet'
