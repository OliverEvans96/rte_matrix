# this script will take an input matrix and solution vector, then solve the system Ax = b
# first, it will compute A = QR for Q unitary and R upper triangular
# next, it will compute inverse(Q)*b which equals R
# finally, R = inv(Q)*B can be solved by back substitution because R is upper triangular

import sys
import time
import numpy as np
import scipy as spy
import scipy.io as scio

# start timing computations
starttime = time.time()
# import A and b from a file given as an argument on the command line when calling this script
mat_contents = scio.loadmat(sys.argv[1])
A = mat_contents['A']
b = mat_contents['b']

# now, the QR algorithm
# first, we compute Q
u1 = A[:,0]
u1 = u1.toarray()
U = np.empty([len(u1),len(u1)])

# set the first column of U as the first column of A
for index in np.arange(0,len(u1),1):
    U[index,0] = u1[index]

# function to project aj onto ui
def proj(ui,aj):
    numerator = np.dot(ui,aj)
    denominator = np.dot(ui,ui)
    return (numerator/denominator)*ui

for column in np.arange(1,len(u1),1):
    projection =  proj(U[:,column-1],A[:,column].toarray())
    for row in np.arange(1,len(u1),1):
        U[row,column] = A[row,column] - projection[row]

Q = np.empty([len(u1),len(u1)])
for column in np.arange(0,len(u1),1):
    Q[:,column] = U[:,column]/np.linalg.norm(U[:,column])

R = np.linalg.inv(Q)@A

X_intermediate = np.linalg.inv(Q)@b
X = np.linalg.inv(R)@X_intermediate
endtime = time.time()
comp_time = endtime - starttime
print(X)
print(comp_time)
