# Hello World

#without fopenmp

```
gcc OMP-hello.c -o OMP-helloOM -fopenmp

./OMP-helloOM

gcc OMP-hello.c -o OMP-hello

./OMP-hello

export OMP_NUM_THREADS=3

./OMP-hello
```

# Execution script

```
#!/bin/bash

#SBATCH -c 8

./exec.sh

```

# Function
```
gcc OMP-hello-function.c -o OMP-hello-function -fopenmp
OMP-hello-function./OMP-hello-function
```

# Thread Affinity
```
export OMP_NUM_THREADS=10
export GOMP_CPU_AFFINITY=0-10
./fib-taskOM

export OMP_PROC_BIND=spread
./fib-taskOM
export OMP_PROC_BIND=close
./fib-taskOM
```

# Data Sharing Clauses
```
OMP-hello-function
gcc OMP-hello-PR-variable.c -o OMP-hello-PR-variable -fopenmp
./OMP-hello-PR-variable
```

# Loop Worksharing
```
gcc OMP-loop-WorkSharing.c -o OMP-loop-WorkSharing -fopenmp
```

# Performance Comparison

Serial version:

```
gcc OMP-matrix-sum.c -o OMP-matrix-sum

time ./OMP-matrix-sum
```

Parallel version:

```
gcc OMP-matrix-sum.c -o OMP-matrix-sumOM -fopenmp

time ./OMP-matrix-sumOM
```

# Race condition

```
gcc OMP-race.c -o OMP-race -fopenmp
./OMP-race
```

# Race condition - synchronization

critical, atomic, ordered

```
gcc OMP-race.c -o OMP-race -fopenmp
./OMP-race
```

# Synchronization

```
gcc OMP-sync.c -o OMP-sync -fopenmp

./OMP-sync

```
# Task

```
gcc OMP-task.c -o OMP-task

./OMP-task 
```

# Fibonacci

```
gcc fib-task.c -o fib-task 

./fib-task 

gcc fib-task.c -o fib-taskOM -fopenmp

./fib-taskOM
```
