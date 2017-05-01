# testing scipy's sparse linear algebra module's factorized() method to solve Ax = b for large sparse matrices

import sys
import time
import numpy as np
import scipy.io as scio
import scipy.sparse.linalg as sla

# start timing computations
start_time = time.time()

# import A and b from file given on command line after script name
mat_contents = scio.loadmat(sys.argv[1])
A = mat_contents['A']
b = mat_contents['b']
b = b.toarray()

# use the factorized method to factorize A and return a function which will solve Ax = b
# (I believe it uses LU factorization)
solve = sla.factorized(A)

X = solve(b)

# stop timing computations
end_time = time.time()
comp_time = end_time - start_time

print(X)
print(comp_time)
