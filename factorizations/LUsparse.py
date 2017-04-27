# uses scipy's sparse linear algebra module to do LU factorization

import sys
import time
import numpy as np
import numpy.linalg as nla
import scipy.io as scio
import scipy.sparse.linalg as sla
import IPython

# start timing computations
start_time = time.time()

# import A and b from file given on command line after script name
mat_contents = scio.loadmat(sys.argv[1])
A = mat_contents['A']
b = mat_contents['b']
b = b.toarray()

# now we simply use the module to solve the equation Ax = b:
Ainv = sla.splu(A)

X = Ainv.solve(b)

# stop timing computations
end_time = time.time()
comp_time = end_time - start_time

print(X)
print(comp_time)
