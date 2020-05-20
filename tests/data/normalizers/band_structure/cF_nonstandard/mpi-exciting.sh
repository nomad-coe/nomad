#!/bin/bash
#
#Job name
#PBS -N Si
#number of nodes and limit of the execution time (hours:minutes:seconds)
#PBS -l nodes=1:ppn=8,walltime=10:00:00
#PBS -W x=NACCESSPOLICY:SINGLEJOB
#queue 
#PBS -q smallRam 

ulimit -s unlimited
export PATH=$PATH:/app/intel/bin/:/app/intel/impi/4.0.3.008/intel64/bin/
export LD_LIBRARY_PATH=/app/intel/lib/intel64/:/app/intel/mkl/lib/intel64/:/app/intel/impi/4.0.3.008/intel64/lib/
cd $PBS_O_WORKDIR

date
export OMP_NUM_THREADS=4
mpirun -n 8 /scratch/dnabok/CODES/EXCITING/carbon-dev/bin/excitingmpi
date
