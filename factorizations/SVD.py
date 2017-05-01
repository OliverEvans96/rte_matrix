# Use SVD Decomposition to solve Ax = b for large sparse matrices

import sys
import time
import numpy as np
import scipy.io as scio

def SVD(matfile):
    # start logging computation time
    start_time = time.time()
    
    # import A and b from the filename specified on the command line when calling SVD.py
    mat_contents = scio.loadmat(matfile)
    A = mat_contents['A']
    b = mat_contents['b']
    A = A.todense()
    
    a1 = A[:,0]
    n = len(a1)
    
    # now, the SVD decomposition
    # first, we compute AA*
    A_Aadj = A@np.transpose(A)
    [w,v] = np.linalg.eig(A_Aadj)
    
    sigma = []
    for ev in range(len(w)):
        sigma.append(np.sqrt(ev))
    sigma.pop(0) # remove the trivial eval zero
    
    # create SIGMA as the diagonal matrix of the singular values of A
    SIGMA = np.diag(sigma)
    
    # create V with the normalized evecs for the evals of A
    V = np.empty([len(a1),len(sigma)])
    vi_norm = v[:,1]/np.linalg.norm(v[:,1])
    
    for column in range(len(sigma)-1):
        for row in range(len(sigma)-1):
            V[row,column] = vi_norm[row]
        vi_norm = v[:,column+2]/np.linalg.norm(v[:,column+2])
    
    
    U = np.empty([len(a1),len(sigma)])
    
    for column in range(len(sigma)-1):
        U[:,column] = (A@V[:,column])/sigma[column]
    
    
    # compute the solution by solving USigmaV*x = b
    
    X1 = np.transpose(U)@b
    X2 = SIGMA@X1
    X = V@X2
    
    # stop logging computation time
    end_time = time.time()
    comp_time = end_time - start_time
    
    X_and_comptime = [X, comp_time]
    
    return(X_and_comptime)

#ans = SVD(sys.argv[1])

#print(ans[0])
#print(ans[1])
