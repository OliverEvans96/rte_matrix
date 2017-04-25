# File Name: block_sor.py
# Description: Try to solve RTE system with block SOR
# Author: Oliver Evans
# Created: Mon Apr 24, 2017 | 03:49pm EDT
# Last Modified: Tue Apr 25, 2017 | 07:33am EDT

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
import scipy.linalg as la
import gen_matrix_2d as gm2
import IPython

# Load system
nx = 50
ny = 50
nth = 16

var_order = [2,1,0]
var_lengths = [nx,ny,nth]

scenario = gm2.KelpScenario()
scenario.load_rte_system_mat('../mat/kelp1_50x50x32_210.mat',var_order)

# Hyperparameters
tol = 1e-3
maxiter = 10000

# The 210 structure has a tightly banded diagonal block,
# and every off-diagonal block is diagonal.
# We'll use block SOR with banded solver after 
# row reduction for diagonal blocks.
# Also, compare banded solver to direct sparse solver.
# Use BLAS DGBMV (banded mat-vec-mult) w/ ku=kl=0
# for off-diagonal multiplication

# The matrix contains nth blocks of size ny, 
# which are themselves blocks of size nx

# Solve
# A[r,r] @ u[r,1] = b[r] - sum_{s!=r} (A[r,s] @ u[s,0])
# Where all of the above indices are blocks except u[,(0/1)],
# which are previous and current iteration respectively

# Block form of Gauss-Seidel iteration
def block_gs(AA,bb,var_lengths,tol,maxiter):
    nx, ny, nth = var_lengths
    nn = AA.shape[0]
    uu = np.zeros([nn,2])

    # Given 1 or 2 theta block indices, return slice of aa
    # in first 1 or 2 dimensions
    def get_block(aa,ind1,ind2=None):
        nn = aa.ndim

        # Length of each block
        bl_len = nx * ny

        # Slice first dimension
        sl1 = slice(bl_len*ind1, bl_len*(ind1+1))

        if ind2==None:
            return aa[sl1]
        else:
            # Slice second dimension
            sl2 = slice(bl_len*ind2, bl_len*(ind2+2))
            return aa[sl1,sl2]

    # Set values in a block
    def set_block(aa,cc,ind1,ind2=None):
        nn = aa.ndim

        # Length of each block
        bl_len = nx * ny

        # Slice first dimension
        sl1 = slice(bl_len*ind1, bl_len*(ind1+1))

        if ind2==None:
            aa[sl1] = cc
            return
        else:
            # Slice second dimension
            sl2 = slice(bl_len*ind2, bl_len*(ind2+2))
            aa[sl1,sl2] = cc
            return

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

    # Loop through iterations
    for it in range(maxiter):

        print("it={:05d}".format(it))

        # Loop through blocks
        for rr in range(nth):
            print("rr={}".format(rr))
            rhs = get_block(bb,rr)
            A_diag = get_block(AA,rr,rr).toarray()

            # Loop over non-diagonal elements to costruct block RHS
            for ss in range(nth):
                # ss != rr
                if ss == rr:
                    continue

                # A[rr,ss] is diagonal, so just multiply accordingly
                # rather than full matrix multiplication
                # Also, use first column of uu[ss] = previous iteration
                rhs_r = (get_block(AA,rr,ss).diagonal()
                        * get_block(uu[:,0],ss))

                # Have to convert to column vector before subtracting
                # Otherwise, 2D row-column broadcast will occur
                # Append dimension
                rhs -= np.expand_dims(rhs_r,-1)

            # Row reduce to eliminate diagonal subblock (ny-1,ny-3)
            # (third from last - one before penultimate)
            # using (ny-2,ny-3)
            # Actually, this wouldn't improve sparsity pattern due to 
            # periodic x bc, so just skip row reduction

            # Solve using LAPACK banded solver
            # If rr < nth/2, then kl = 2*nx
            # otherwise ku = 2*nx
            # other bandwidth is just nx


            if rr < nth / 2:
                kl = 2*nx
                ku = nx
            else:
                kl = nx
                ku = 2*nx

            # First, need to create aa_band in proper banded format
            aa_band = arr2band(A_diag,kl,ku)

            # Use LAPACK double general banded solver
            lub, piv, xx, info = la.lapack.dgbsv(kl,ku,aa_band,rhs)

            # Save result to rr block of second column of u
            # xx is (nn x 1), but it seems necessary to convert
            # it to 1d before using it to set uu
            set_block(uu[:,1],xx[:,0],rr)

        # Calculate error (difference from last iteration)
        err = np.sum(abs(uu[:,0]-uu[:,1])) / nn
        print("err = {:.3e}".format(err))

        # Terminate iteration if error criteria is met
        if err < tol:
            return uu[:,1]

        # Update to next iteration
        uu[:,0] = uu[:,1]

block_gs(scenario._rte_matrix,scenario._rte_rhs,var_lengths,tol,maxiter)

