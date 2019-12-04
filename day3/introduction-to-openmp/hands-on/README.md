# Matrix Transposition

The following code was used in the book "lots of core"

Parallelize the Transposition function at "hands-on/Transposition/Transpose.cc"

compilation:
```
cd hands-on/Transposition
make clean 
make runme-CPU
```

Execution:
```
./runme-CPU [size] [repetition]
```

Parameter [size] is size of array
Parameter [trials] define how many times the operations will be performed

Example: ./runme-CPU 3000 100

# IronBar 

Parallelize this code using OpenMP.

This code was used as warmup for Marathon of Parallel Programming at ERAD-SP 2018 (http://lspd.mackenzie.br/marathon/18/warmup.pdf)

compilation:
```
gcc -O3 -fopenmp ironbar.c -o ironbar
```
Execution:
```
./ironbar < input2
```

# Harmonic Progression SUM (WSCAD 2016 warmup)​

Parallelize this code using OpenMP.

This code implements the Harmonic Progression sum and was used as warmup for Marathon of Parallel Programming at WSCAD 2016 (http://www.wscad-2016.ufs.br)

Compilation:
```
cd sum
g++ sum.cpp -o sum -fopenmp
```

Execution:

Fast example:
```
./sum < sum.in
```

Slower example:
```
./sum < sum2.in
```

# QuickSort algorithm

Paralelize this code using OpenMP Tasks.

Compilation:
```
gcc -fopenmp -o quicksort quicksort.c
```
	
Execution:
```
    ./quicksort [size] 
```
Parameter [size] is size of array

Example: ./quicksort 20000000

# Nbody

The following code was used in the book "lots of core"

this code simulate particle interaction according to Newton's Law

compilation:
```
g++ -O3 -fopenmp -ffast-math -march=native -mtune=native -o nbody-v0s nbody-v0.cc
```

execution:
```
./nbody-v0s
```

In order to change the size of problem go to line (101):
```
    // Number of particles (bodies) in our simulated "universe"
    const size_t n = 50000;
```
