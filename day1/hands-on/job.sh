#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=8 
#SBATCH --time=00:30:00
##SBATCH --reservation=treinamento_2 

#module load intel/mpi
#mpirun /bin/hostname

#ls -ltr 
/bin/hostname
./prog.x
sleep 100
