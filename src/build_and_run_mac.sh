#!/bin/sh  

rm ./les_main_mpi
procPerRow=2
procPerCol=2

scons -f SConstruct.mac wv_debug=1 ocl=0 mpi=1 procPerRow=${procPerRow} procPerCol=${procPerCol}

mpiexec-openmpi-gcc49 -np $((procPerRow*procPerCol)) ./les_main_mpi 
