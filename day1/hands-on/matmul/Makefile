#simple makefile for the matmul program
FC = gcc
FFLAGS = -O3 -I/opt/intel/composer_xe_2013_sp1.2.144/mkl/include

LIB = -fopenmp -Wl,--rpath,/opt/intel/composer_xe_2013_sp1.2.144/mkl/lib/intel64/ -L/opt/intel/composer_xe_2013_sp1.2.144/mkl/lib/intel64/ -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm

LIB2 = -fopenmp -Wl,--rpath,/opt/intel/composer_xe_2013_sp1.2.144/mkl/lib/intel64/ -L/opt/intel/composer_xe_2013_sp1.2.144/mkl/lib/intel64/ -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lm

default: matmul matmul_threads

matmul : matmul.c
	$(FC) $(FFLAGS) matmul.c -o matmul.x $(LIB)

matmul_threads : matmul.c
	$(FC) $(FFLAGS) matmul.c -o matmul_threads.x $(LIB2)

clean: 
	- \rm -f *.x  


