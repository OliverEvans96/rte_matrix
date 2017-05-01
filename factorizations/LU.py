# use LU decomp to solve system Ax = b for matrices in rte_matrix/

import sys
import time
import numpy as np
import scipy.io as scio

def LU(matfile):

    # start logging computation time
    start_time = time.time()
    
    # import A and b from command line argument after script name
    mat_contents = scio.loadmat(matfile)
    A = mat_contents['A']
    b = mat_contents['b']
    
    # now, the LU decomposition
    # create L_n (the matrix that eliminates entries below the main diag of A in the nth row)
    A = A.todense()
    a1 = A[:,0]
    
    L = np.empty([len(a1),len(a1)])
    for column in range(len(a1)):
        for row in np.arange(column,len(a1),1):
            if row == column:
                L[row,column] = 1
            else:
                L[row,column] = -1*A[row,column]/A[column,column]
    
    U = L@A
    
    X1 = np.linalg.inv(L)@b
    X = np.linalg.inv(U)@X1

    # stop logging computation time
    end_time = time.time()
    comp_time = end_time - start_time

    X_and_comptime = [X, comp_time]
    return(X_and_comptime)

# to run this script on its own, un-comment the following lines
#ans = LU(sys.argv[1])
#print(ans[0])
#print(ans[1])
