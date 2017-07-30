#!/bin/sh  

rm ./les_main

scons -f SConstruct.mac v=0 gr_debug=0 wv_debug=0 mpi_new_wv=0 nested=0 ocl=0 mpi=0 procPerRow=1 procPerCol=1 $*

time ./les_main
