# File Name: block_sor.py
# Description: Try to solve RTE system with block SOR
# Author: Oliver Evans
# Created: Mon Apr 24, 2017 | 03:49pm EDT
# Last Modified: Mon May 01, 2017 | 10:22pm EDT

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
import matplotlib.pyplot as plt

plt.ion()

# Grid parameters
nx = 20
ny = 20
nth = 24
nrows = nx * ny * nth

var_order = [2,1,0]
var_lengths = [nx,ny,nth]

# IOPs
def surf_bc_fun(th):
    return np.sin(th)

# Normalized mock VSF
def vsf(th):
    return 2 * np.exp(-th/2) / (1 - np.exp(-np.pi/2))

abs_water = 1
sct_water = 1
abs_kelp = 5
sct_kelp = 1
iops = [vsf,abs_water,sct_water,abs_kelp,sct_kelp]

# Load matrix to test
scenario = gm2.KelpScenario(surf_bc_fun,iops)
scenario.set_num_grid_points(nx,ny,nth)
scenario.load_rte_system_mat('../mat/kelp1_{}x{}x{}_210.mat'
        .format(nx,ny,nth),var_order)

dx = 1/nx
dy = 1/ny

# Maximum number density of individuals 
ind_up = 2

# Max/min abs. coef.
a_up = abs_kelp * ind_up
a_dn = abs_water
# Max/min scat. coef.
b_up = sct_kelp * ind_up
b_dn = sct_water

# Hyperparameters
tol = 1e-8
maxiter = 500

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
def block_gs(AA,bb,var_lengths,tol,maxiter,x0=None,true_sol=None):
    nx, ny, nth = var_lengths
    nn = AA.shape[0]

    # Initialize error
    err = np.zeros(maxiter)

    uu = np.zeros([nn,2])
    # Initial guess f provided
    if x0 is not None:
        uu[:,0] = x0

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
            sl2 = slice(bl_len*ind2, bl_len*(ind2+1))
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

    # At each iteration, make the matrix diagonally dominant
    # by adding mu*u[:,1] to LHS and mu*u[:,0] to RHS.
    # After many iterations, the two should be approx. equal.
    mu = 1/dx + 4/dy + b_up * (2*np.pi - 1) - a_dn

    # Loop through iterations
    for it in range(maxiter):

        print("it={:05d}".format(it))

        # Loop through blocks
        for rr in range(nth):
            # Get diagonal block
            #print("rr={}".format(rr))
            rhs = get_block(bb,rr).toarray()
            A_diag = get_block(AA,rr,rr).toarray()

            # Create diagonal dominance
            A_diag += mu * np.eye(nx*ny)

            # Adjust RHS accordingly
            rhs += mu * np.expand_dims(get_block(uu[:,0],rr),-1)

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
                # Otherwise, 2D row-column broadcast will occur.

                # Append dimension
                rhs -= np.expand_dims(rhs_r,-1)

            # Solve using LAPACK banded solver
            # If rr < nth/2, then kl = 2*nx
            # otherwise ku = 2*nx
            # other (ku/kl) bandwidth is just nx

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
        if true_sol is None:
            err[it] = np.sum(abs(uu[:,0]-uu[:,1])) / nn
        else:
            err[it] = np.sum(abs(true_sol-uu[:,1])) / nn
        print("err = {:.3e}".format(err[it]))
        print()

        # Terminate iteration if error criteria is met
        if err[it] < tol:
            return (uu[:,1],err[:it+1])

        # Update to next iteration
        uu[:,0] = uu[:,1]

    print("Failed to converge to {} after {} iterations".format(tol,maxiter))
    return (uu[:,1],err)

# Plot direct solution
#scenario._rte_sol = block_gs(scenario._rte_matrix,scenario._rte_rhs,var_lengths,tol,maxiter)
scenario.solve_system()
scenario.reshape_rad()
scenario.calc_irrad()
scenario.plot_irrad('tmp.png')
true_sol = scenario._rte_sol
gm2.plt.show()

# Plot iterative solution every np iterations up to maxiter
maxiter = 200
nstep = 50
nplots = maxiter // nstep
err = np.zeros(maxiter)

x1 = np.zeros(nx*ny*nth)
for it in range(nplots):
    scenario._rte_sol,tmperr = block_gs(scenario._rte_matrix,scenario._rte_rhs,var_lengths,tol,nstep,x0=x1,true_sol=true_sol)

    # Update initial guess
    x1 = scenario._rte_sol

    # Append error
    err[it*nstep:(it+1)*nstep] = tmperr

    plt.figure(2*it+1)
    # Plot error 
    plt.plot(err)
    plt.title('error')

    # Plot 
    scenario.reshape_rad()
    scenario.calc_irrad()
    plt.figure(2*it+2)
    scenario.plot_irrad('tmp.png')

    plt.show()

