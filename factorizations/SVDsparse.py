# uses scipy's sparse linear algebra module to do SVD factorization to solve equations Ax = b

import sys
import time
import numpy as np
import scipy.io as scio
import scipy.sparse.linalg as sla

def SVDsparse(matfile):
    # start timing computations
    start_time = time.time()
    
    # import A and b from file given on command line after script name
    mat_contents = scio.loadmat(matfile)
    A = mat_contents['A']
    b = mat_contents['b']
    b = b.toarray()
    
    # compute singular values and their corresponding vectors using scipy.sparse.linalg
    [U,sigma,Vt] = sla.svds(A,1599)
    SIGMA = np.diag(sigma)
    
    X1 = np.transpose(U)@b
    X2 = SIGMA@X1
    X = np.transpose(Vt)@X2
    
    # stop timing computations
    end_time = time.time()
    comp_time = end_time - start_time
    
    X_and_comptime = [X, comp_time]
    
    return(X_and_comptime)


#ans = SVDsparse(sys.argv[1])

#print(ans[0])
#print(ans[1])
