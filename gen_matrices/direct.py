# File Name: band.py
# Description: Direct methods for linear system
# Created: Thu May 04, 2017 | 01:54am EDT
# Last Modified: Thu May 04, 2017 | 07:17am EDT

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#                           GNU GPL LICENSE                            #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#                                                                      #
# Copyright Oliver Evans 2017 <oliverevans96@gmail.com>                #
#                                                                      #
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with this program. If not, see <http://www.gnu.org/licenses/>. #
#                                                                      #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

import numpy as np
import scipy.sparse.linalg as sla
import time

# Banded
# spsolve
# LUsparse
# SVDsparse

# numpy array -> LAPACK banded
# Don't try to convert straight from csr. It's very slow at slicing.
# http://www.netlib.org/lapack/explore-html-3.4.2/d3/d49/group__double_g_bsolve.html#gafa35ce1d7865b80563bbed6317050ad7
# http://www.netlib.org/lapack/lug/node124.html
def arr2band(aa_arr,kl,ku):
    nn = aa_arr.shape[0]
    ldab = 2 * kl + ku + 1
    aa_band = np.zeros([ldab,nn])
    
    for jj in range(nn):
        #print("jj/nn={}/{}".format(jj,nn))
        for ii in range(max(0,jj-ku),min(nn,jj+kl)):
            aa_band[kl+ku+ii-jj,jj] = aa_arr[ii,jj]

    return aa_band

def band_slv(A,b,kl,ku):
    t1 = time.time()
    sol = la.solve_banded(arr2band(A,kl,ku),b)
    t2 = time.time()

    return (sol, t1, t2)

def spsolve(A,b):
    t1 = time.time()
    sol = sla.spsolve(A,b)
    t2 = time.time()

    return (sol, t1, t2)

def LUsparse(A, b):
    # start timing computations
    start_time = time.time()
    
    # now we simply use the module to solve the equation Ax = b:
    Ainv = sla.splu(A)
    
    X = Ainv.solve(b)
    
    # stop timing computations
    end_time = time.time()
    comp_time = end_time - start_time

    return (sol, start_time, end_time)

def SVDsparse(A, b):
    # start timing computations
    start_time = time.time()
    
    # compute singular values and their corresponding vectors using scipy.sparse.linalg
    [U,sigma,Vt] = sla.svds(A,A.shape[0])
    SIGMA = np.diag(sigma)
    
    X1 = np.transpose(U)@b
    X2 = SIGMA@X1
    X = np.transpose(Vt)@X2
    
    # stop timing computations
    end_time = time.time()
    comp_time = end_time - start_time
    
    X_and_comptime = [X, comp_time]
    
    return(X_and_comptime)

# Collect all direct methods into a dictionary
direct_method = dict(
    band = band_slv,
    lu = LUsparse,
    svd = SVDsparse)

