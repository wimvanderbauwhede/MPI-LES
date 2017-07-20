#!/bin/sh  

rm ./les_main_mpi
procPerRow=4
procPerCol=4

scons -f SConstruct.mac v=0 wv_debug=0 nested=1 ocl=0 mpi=1 procPerRow=${procPerRow} procPerCol=${procPerCol}

mpiexec-openmpi-gcc49 -np $((procPerRow*procPerCol)) ./les_main_mpi 
