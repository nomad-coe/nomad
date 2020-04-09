#!/bin/bash
#
#This is an example script for running exciting 
#
#Job name
#PBS -N Si-RMTbig-0.9795861087156
#number of nodes and limit of the execution time (hours:minutes:seconds)
#PBS -l nodes=1:ppn=32,walltime=100:00:00
#queue 
#PBS -q smallRam

cd $PBS_O_WORKDIR
export PATH=$PATH:/app/intel/bin/:/app/intel/impi/4.0.3.008/intel64/bin/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/app/intel/lib/intel64/:/app/intel/mkl/lib/intel64/
export OMP_NUM_THREADS=2

# mpirun -n 16 <program> launches a parallel program with 16 processes
# pay attention how the number of processes compares with nodes*ppn above
ulimit -s unlimited 

hostname
date
mpirun -n 16 -perhost 16 /scratch/gulans/mortadella/latest/exciting/bin/excitingmpi 
date
