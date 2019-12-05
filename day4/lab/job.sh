#!/bin/bash
#SBATCH -N 2
#SBATCH --ntasks-per-node=8 
#SBATCH --time=00:30:00
##SBATCH --reservation=treinamento_2 

module load intel/mpi
srun ./mpi_hello.x
