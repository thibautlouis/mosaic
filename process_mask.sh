#!/bin/sh
export LIBRARY_PATH="$HOME/local/lib:/opt/intel/compilers_and_libraries_2017.4.181/mac/compiler/lib:$MKLPATH:/opt/local/lib/"
export LD_LIBRARY_PATH="$LIBRARY_PATH"

process_mask distance.dict
