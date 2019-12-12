################
# Compiler settings and flags
CC=gcc
OPT=-O3 -march=native -ffast-math -fexpensive-optimizations -ftree-vectorize -Wall
EXE= ljmd-noopt.x ljmd-orig.x ljmd-opt.x \
	ljmd-omp-atomic.x ljmd-omp-reduce.x \
	ljmd-mpi-atomic.x ljmd-mpi-reduce.x \
	ljmd-cell-atomic.x ljmd-cell-reduce.x \
	ljmd-full-atomic.x ljmd-full-reduce.x

################

default: $(EXE)

ljmd-noopt.x: ljmd-orig.c
	$(CC) -o $@ -O0 -g $< -lm

ljmd-orig.x: ljmd-orig.c
	$(CC) -o $@ $(OPT) -g $< -lm

ljmd-opt.x: ljmd-opt.c
	$(CC) -o $@ $(OPT) -g $< -lm

ljmd-omp-atomic.x: ljmd-omp-atomic.c
	$(CC) -o $@ -fopenmp $(OPT) -g $< -lm

ljmd-omp-reduce.x: ljmd-omp-reduce.c
	$(CC) -o $@ -fopenmp $(OPT) -g $< -lm

ljmd-mpi-atomic.x: ljmd-mpi-atomic.c
	mpicc -o $@ -fopenmp $(OPT) -g $< -lm

ljmd-mpi-reduce.x: ljmd-mpi-reduce.c
	mpicc -o $@ -fopenmp $(OPT) -g $< -lm

ljmd-cell-atomic.x: ljmd-cell-atomic.c
	$(CC) -o $@ -fopenmp $(OPT) -g $< -lm

ljmd-cell-reduce.x: ljmd-cell-reduce.c
	$(CC) -o $@ -fopenmp $(OPT) -g $< -lm

ljmd-full-atomic.x: ljmd-full-atomic.c
	mpicc -o $@ -fopenmp $(OPT) -g $< -lm

ljmd-full-reduce.x: ljmd-full-reduce.c
	mpicc -o $@ -fopenmp $(OPT) -g $< -lm


clean:
	rm -f *.x *.dat *.xyz *.o *.mod perf.data *~

bench: bench1 bench2 bench3

bench1: $(EXE)
	@for s in noopt orig opt; do \
	    echo "./ljmd-$$s.x < argon_108.inp"; \
	    time ./ljmd-$$s.x < argon_108.inp >> /dev/null; \
	done
	@for t in 1 2 4; \
	    do for s in omp-atomic omp-reduce cell-atomic cell-reduce; do \
	    echo "env OMP_NUM_THREADS=$$t ./ljmd-$$s.x < argon_108.inp"; \
	    time env OMP_NUM_THREADS=$$t ./ljmd-$$s.x < argon_108.inp >> /dev/null; \
	done; done
	@for t in 1 2 4; \
	    do for s in omp-atomic omp-reduce cell-atomic cell-reduce; do \
	    echo "env OMP_NUM_THREADS=$$t ./ljmd-$$s.x < argon_108.inp" ; \
	    time env OMP_NUM_THREADS=$$t ./ljmd-$$s.x < argon_108.inp >> /dev/null; \
	done; done
	@m=1; for t in 1 2 4; \
	    do for s in mpi-atomic mpi-reduce full-atomic full-reduce; do \
	    echo "mpirun -np $$m env OMP_NUM_THREADS=$$t ./ljmd-$$s.x < argon_108.inp"; \
	    time mpirun -np $$m env OMP_NUM_THREADS=$$t ./ljmd-$$s.x < argon_108.inp >> /dev/null; \
	done; done
	@m=2; for t in 1 2; \
	    do for s in mpi-atomic mpi-reduce full-atomic full-reduce; do \
	    echo "mpirun -np $$m env OMP_NUM_THREADS=$$t ./ljmd-$$s.x < argon_108.inp"; \
	    time mpirun -np $$m env OMP_NUM_THREADS=$$t ./ljmd-$$s.x < argon_108.inp >> /dev/null; \
	done; done
	@m=4; for t in 1; \
	    do for s in mpi-atomic mpi-reduce full-atomic full-reduce; do \
	    echo "mpirun -np $$m env OMP_NUM_THREADS=$$t ./ljmd-$$s.x < argon_108.inp"; \
	    time mpirun -np $$m env OMP_NUM_THREADS=$$t ./ljmd-$$s.x < argon_108.inp >> /dev/null; \
	done; done

bench2:
	@for s in noopt orig opt; do \
	    echo "./ljmd-$$s.x < argon_2916.inp"; \
	    time ./ljmd-$$s.x < argon_2916.inp >> /dev/null; \
	done
	@for t in 1 2 4; \
	    do for s in omp-atomic omp-reduce cell-atomic cell-reduce; do \
	    echo "env OMP_NUM_THREADS=$$t ./ljmd-$$s.x < argon_2916.inp"; \
	    time env OMP_NUM_THREADS=$$t ./ljmd-$$s.x < argon_2916.inp >> /dev/null; \
	done; done
	@for t in 1 2 4; \
	    do for s in omp-atomic omp-reduce cell-atomic cell-reduce; do \
	    echo "env OMP_NUM_THREADS=$$t ./ljmd-$$s.x < argon_2916.inp" ; \
	    time env OMP_NUM_THREADS=$$t ./ljmd-$$s.x < argon_2916.inp >> /dev/null; \
	done; done
	@m=1; for t in 1 2 4; \
	    do for s in mpi-atomic mpi-reduce full-atomic full-reduce; do \
	    echo "mpirun -np $$m env OMP_NUM_THREADS=$$t ./ljmd-$$s.x < argon_2916.inp"; \
	    time mpirun -np $$m env OMP_NUM_THREADS=$$t ./ljmd-$$s.x < argon_2916.inp >> /dev/null; \
	done; done
	@m=2; for t in 1 2; \
	    do for s in mpi-atomic mpi-reduce full-atomic full-reduce; do \
	    echo "mpirun -np $$m env OMP_NUM_THREADS=$$t ./ljmd-$$s.x < argon_2916.inp"; \
	    time mpirun -np $$m env OMP_NUM_THREADS=$$t ./ljmd-$$s.x < argon_2916.inp >> /dev/null; \
	done; done
	@m=4; for t in 1; \
	    do for s in mpi-atomic mpi-reduce full-atomic full-reduce; do \
	    echo "mpirun -np $$m env OMP_NUM_THREADS=$$t ./ljmd-$$s.x < argon_2916.inp"; \
	    time mpirun -np $$m env OMP_NUM_THREADS=$$t ./ljmd-$$s.x < argon_2916.inp >> /dev/null; \
	done; done

bench3:
	@for s in noopt orig opt; do \
	    echo "./ljmd-$$s.x < argon_78732.inp"; \
	    time ./ljmd-$$s.x < argon_78732.inp >> /dev/null; \
	done
	@for t in 1 2 4; \
	    do for s in omp-atomic omp-reduce cell-atomic cell-reduce; do \
	    echo "env OMP_NUM_THREADS=$$t ./ljmd-$$s.x < argon_78732.inp"; \
	    time env OMP_NUM_THREADS=$$t ./ljmd-$$s.x < argon_78732.inp >> /dev/null; \
	done; done
	@for t in 1 2 4; \
	    do for s in omp-atomic omp-reduce cell-atomic cell-reduce; do \
	    echo "env OMP_NUM_THREADS=$$t ./ljmd-$$s.x < argon_78732.inp" ; \
	    time env OMP_NUM_THREADS=$$t ./ljmd-$$s.x < argon_78732.inp >> /dev/null; \
	done; done
	@m=1; for t in 1 2 4; \
	    do for s in mpi-atomic mpi-reduce full-atomic full-reduce; do \
	    echo "mpirun -np $$m env OMP_NUM_THREADS=$$t ./ljmd-$$s.x < argon_78732.inp"; \
	    time mpirun -np $$m env OMP_NUM_THREADS=$$t ./ljmd-$$s.x < argon_78732.inp >> /dev/null; \
	done; done
	@m=2; for t in 1 2; \
	    do for s in mpi-atomic mpi-reduce full-atomic full-reduce; do \
	    echo "mpirun -np $$m env OMP_NUM_THREADS=$$t ./ljmd-$$s.x < argon_78732.inp"; \
	    time mpirun -np $$m env OMP_NUM_THREADS=$$t ./ljmd-$$s.x < argon_78732.inp >> /dev/null; \
	done; done
	@m=4; for t in 1; \
	    do for s in mpi-atomic mpi-reduce full-atomic full-reduce; do \
	    echo "mpirun -np $$m env OMP_NUM_THREADS=$$t ./ljmd-$$s.x < argon_78732.inp"; \
	    time mpirun -np $$m env OMP_NUM_THREADS=$$t ./ljmd-$$s.x < argon_78732.inp >> /dev/null; \
	done; done

.PHONY: clean bench bench1 bench2 bench3
