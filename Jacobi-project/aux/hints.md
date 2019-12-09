## MPI - 1D Decomposition

### Assignments

The parameters of the algorithm are such:


1.  The grid matrix must be completely distributed, no replicating the
    matrix on all processors. In this exercise, only use a 1 dimensional decomposition (see
    [Figure 2](#Figure_2)).
    
    ![Figure 2](jacobiFigure2.jpg)
    
2.  The whole process must be parallel, that includes initialization of
    the grid and boundary conditions, the iterative evolution and the final dump on file of the resulting grid. 
    
3.  Implement an efficient data exchange between processes.
    
4.  Handle dimensions even if not multiple of the number of processes.

Here is a guideline for the process that parallel programmers use to do
this:

1.  Study the serial algorithm and see where parallelism can be
    exploited. Also think about how the data can be divided. Best way
    to do this is on a piece of paper, drawing out the layout
    conceptually before you even touch the code.
    
2.  Still on paper, figure out how this conceptualization moves to being
    expressed in the parallel programming language you want to use.
    What MPI calls do you need to use? Which processors will be doing
    what work? STILL ON PAPER.
    
3.  Now begin programming the algorithm up in MPI.
    
4.  Test the program on a small matrix and processor count to make sure
    it is doing what you expect it to do.
    
5.  Once you are satisfied it works, scale it up.

With this in mind, go through this process to implement a 1-D
decomposition of the Jacobi iteration algorithm.


### Tips

-   To set up the initial matrix, you will need to figure out which
    values go in what chunk of the distributed matrix. Think carefully
    about the data that each parallel chunk of work needs to have on
    it.
    
-   Notice the value of each matrix element depends on the adjacent
    elements from the previous matrix. In the distributed matrix, this
    has consequences for the boundary elements, in that if you
    straightforwardly divide the matrix up by rows, elements that are
    needed to compute a matrix element will reside on a different
    processor. Think carefully about how to allocate the piece of the
    matrix on the current processor, and what communication needs to be
    performed before computing the matrix elements. [Figure
    2](#Figure_2). is an illustration of one communication patter that
    can be used.
    
    
-   It is requested to write a function that will print the
    distributed matrix, so that you have the ability to check to see
    if things are going the way you want.

-   To perform a data exchange with a “dummy” process you can use
    [MPI_PROC_NULL](http://mpi-forum.org/docs/mpi-1.1/mpi-11-html/node53.html)

-   A reference of MPI routines can be found at:
    <http://mpi-forum.org/docs/mpi-1.1/mpi-11-html/node182.html>






