#!/usr/bin/env python
import time
from mpi4py import MPI
import iso_dict
import sys
import os

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

print rank, size

p = iso_dict.flipperDict()
p.read_from_file(sys.argv[1])

iStart = p['iStart']
iStop = p['iStop']

delta = (iStop - iStart)/size

if delta == 0:
    raise ValueError, 'Too many processors for too small a  loop!'

print delta


for i in range(iStart+rank, iStop, size):
    print "compiling  in iteration %03d"%i
    os.system('iso_generate_spectra.py %s %03d'%(sys.argv[1],i))



