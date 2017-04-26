# use LU decomp to solve system Ax = b for matrices in rte_matrix/

import sys
import time
import numpy as np
import scipy.io as scio

# start logging computation time
start_time = time.time()

# import A and b from command line argument after script name
mat_contents = scio.loadmat(sys.argv[1])
A = mat_contents['A']
b = mat_contents['b']

