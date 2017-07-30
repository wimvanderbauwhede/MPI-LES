#!/bin/sh  

rm ./les_main

procPerRow=1
procPerCol=1

scons -f SConstruct.mac v=0 gr_debug=0 wv_debug=0 mpi_new_wv=1 nested=0 ocl=0 mpi=0 procPerRow=${procPerRow} procPerCol=${procPerCol} $3

#time ./les_main
