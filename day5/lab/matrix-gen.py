#!/usr/bin/env python3

from numpy import linalg
from numpy import matlib
from numpy import random
from numpy import matmul
from numpy import matrix

random.seed(1543)

amat = matlib.rand(10,10)
bmat = linalg.inv(amat)
dmat = matrix.transpose(bmat)
cmat = matmul(amat,bmat)

def printmat(mat):
    rows = mat.shape[0]
    columns = mat.shape[1]
    print(rows,columns)
    for i in range(rows):
        for j in range(columns):
            print(mat.item((i,j)), end=" ")
        print()

printmat(amat)
printmat(dmat)
