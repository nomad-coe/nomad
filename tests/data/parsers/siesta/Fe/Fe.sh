#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=30:00
#SBATCH --ntasks=2
#SBATCH -p debug
#SBATCH --output=out.log
#SBATCH --error=err.out.log
#SBATCH --exclusive

#load the environment

ulimit -s unlimited
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so
export I_MPI_FABRICS=shm:dapl
module load   intel/2016a    

srun  --resv-ports  -N 1 -n 2 /home/andrea/siesta-4.0/Obj/siesta < Fe.fdf > out