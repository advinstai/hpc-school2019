################
# Compiler settings and flags
FC=gfortran
#FC=mpifort
#OPENMP=-fopenmp
OPENMP=
#FFLAGS= -O0 -Wall -g $(OPENMP)
FFLAGS=-O3 -march=native -ffast-math -fexpensive-optimizations -ftree-vectorize -g $(OPENMP)
################

ljmd.x: ljmd.f90
	$(FC) -o $@ $(FFLAGS) $<

clean:
	rm -f ljmd.x *.dat *.xyz *.o *.mod perf.data

bench: bench1 bench2 bench3

bench1: ljmd.x
	time ./ljmd.x < argon_108.inp

bench2:
	time ./ljmd.x < argon_2916.inp

bench3:
	time ./ljmd.x < argon_78732.inp

.PHONY: clean bench bench1 bench2 bench3
