#!/bin/sh  

rm ./les_main_mpi

procPerRow=$1
procPerCol=$2

scons -f SConstruct.mac v=0 gr_debug=0 wv_debug=0 mpi_new_wv=1 nested=1 ocl=0 mpi=1 procPerRow=${procPerRow} procPerCol=${procPerCol} $3

mpiexec-openmpi-gcc49 -np $((procPerRow*procPerCol)) ./les_main_mpi 
