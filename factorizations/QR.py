# this script will take an input matrix and solution vector, then solve the system Ax = b
# first, it will compute A = QR for Q unitary and R upper triangular
# next, it will compute inverse(Q)*b which equals R
# finally, R = inv(Q)*B can be solved by back substitution because R is upper triangular

import sys
import time
import numpy as np
import numpy.linalg as npla
import scipy.io as scio

def QR(matfile):
    # start timing computations
    start_time = time.time()
    
    # import A and b from a file given as an argument on the command line when calling this script
    mat_contents = scio.loadmat(matfile)
    A = mat_contents['A']
    b = mat_contents['b']
    
    # now, the QR decomposition
    # first, we compute Q
    u1 = A[:,0]
    u1 = u1.toarray()
    U = np.empty([len(u1),len(u1)])
    
    # set the first column of U as the first column of A
    for index in range(len(u1)):
        U[index,0] = u1[index]
    
    # function to project aj onto ui
    def proj(ui,aj):
        numerator = np.dot(ui,aj)
        denominator = np.dot(ui,ui)
        return (numerator/denominator)*ui
    
    for column in range(1,len(u1)):
        projection =  proj(U[:,column-1],A[:,column].toarray())
        U[:,column] = A[:,column].toarray()[:,0] - projection
    
    Q = np.empty([len(u1),len(u1)])
    for column in range(len(u1)):
        Q[:,column] = U[:,column]/npla.norm(U[:,column])
    
    R = npla.inv(Q)@A
    
    X_intermediate = npla.inv(Q)@b
    X = npla.inv(R)@X_intermediate
    
    # stop timing computations
    end_time = time.time()
    comp_time = end_time - start_time
    
    X_and_comptime = [X, comp_time]
    
    return(X_and_comptime)


#ans = QR(sys.argv[1])

#print(ans[0])
#print(ans[1])
