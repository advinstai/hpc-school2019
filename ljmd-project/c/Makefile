################
# Compiler settings and flags
CC=gcc
#CC=mpicc
#OPENMP=-fopenmp
OPENMP=
#CFLAGS= -O0 -Wall -g $(OPENMP)
CFLAGS=-O3 -march=native -ffast-math -fexpensive-optimizations -ftree-vectorize -g $(OPENMP)
################

ljmd.x: ljmd.c
	$(CC) -o $@ $(CFLAGS) $< -lm

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
