# File Name: block_sor.py
# Description: Try to solve RTE system with block SOR
# Author: Oliver Evans
# Created: Mon Apr 24, 2017 | 03:49pm EDT
# Last Modified: Mon Apr 24, 2017 | 09:25pm EDT

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

# Load system
var_order = [2,1,0]
scenario = gm2.KelpScenario(None,None)
scenario.load_rte_system_mat('../mat/kelp1_50x50x32_210.mat',var_order)

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
def block_gs(AA,bb,var_lengths):
    nx, ny, nth = var_lengths
    nn = AA.shape[0]
    uu = np.zeros([nn,2])

    # Given mm theta block indices, return slice of aa
    # in first mm dimensions (need nn <= mm)
    def get_block(aa,*ind_list):
        nn = aa.ndim
        mm = len(ind_list)

        if mm > nn:
            raise IndexError("{} slices given for {}-d array"
                    .format(mm,mm))

        # Length of each block
        bl_len = nx * ny

        # Slice each dimension
        slice_list = [slice(bl_len*bi, bl_len*(bi+1)) for bi in ind_list]

        return aa[sl_list]
    
    for rr in range(nth):
        rhs = get_block(bb,rr)

        for ss in range(nth):
            # ss != rr
            if ss == rr:
                continue

            # A[rr,ss] is diagonal, so just multiply accordingly
            # rather than full, expensive matrix multiplication
            rhs -= diag(get_block(A,rr,ss)) * get_block(u,ss)[:,0]
        






