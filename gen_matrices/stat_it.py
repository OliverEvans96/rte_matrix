# File Name: stat_it.py
# Description: Stationary iterative methods
# Author: Oliver Evans
# Created: Wed May 03, 2017 | 10:14am EDT
# Last Modified: Wed May 03, 2017 | 04:28pm EDT

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

import IPython
import numpy as np
import numpy.linalg as la
from scipy import io
import scipy.sparse as sp
import scipy.sparse.linalg as sla
import matplotlib
# Non-interactive matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time
import gershgorin as gg

def st_itr_tmpl(AA, bb, N_fun,
        tol=1e-4, maxiter=1000,
        true_sol=None, 
        name=None,
        method=None,
        plot_disks=False,
        diskimg=None,
        *args, **kwargs):
    """
    Use notation from 
    https://etd.ohiolink.edu/rws_etd/document/get/akron1342456200/inline

    Stationary iterative method template
    Solve Ax = b
    Let L, D, U be lower, diagonal, & upper triangular 
    Split A = N - L,
    where NN = N_fun(DD,LL,UU)
    true solution, x can be used for error checking if provided
    Otherwise, diff. from prev. iter is used

    If plotdisks=True, M is calculated and Gershgorin disks are
    plotted, and saved in diskimg if specified, otherwise to
    ../img/disks/name_method.png

    *args & **kwargs passed to N_fun

    Return solution, error array, and status
    status = 0 for successful convergence
    status = 1 for maxiter reached
    """

    # Determine matrix leading dimension
    nn = AA.shape[0]

    if sp.issparse(AA):
        # scipy.sparse
        pkg = sp
        # Create diagonal matrix from vector of entries
        diag = sp.diags
        # Which algorithm to use to solve a triangular system
        tri_solve = sla.spsolve
    else:
        # numpy
        pkg = np
        diag = np.diag
        tri_solve = lp.trtrs
        
    # Isolate diagonal & strictly lower & upper triangular parts
    DD = diag(AA.diagonal())
    LL = pkg.tril(AA) - DD
    UU = pkg.triu(AA) - DD

    # Calculate NN, PP
    NN = N_fun(DD,LL,UU,*args,**kwargs)
    PP = NN - AA

    # Convert to column format
    if sp.issparse(AA):
        NN = NN.tocsc()
        PP = PP.tocsc()

    # MM = inv(NN) * PP
    # Only necessary to compute explicity for convergence checking
    # Method converges if & only if sectral radius of M < 1
    if plot_disks:
        # Determine name by matrix size if not specified
        if name is None:
            name = "{}x{}".format(nn,nn)

        print("Calculating iteration matrix")
        t = time.time()
        MM = tri_solve(NN,PP)
        print("(took {:.3e} seconds)".format(time.time()-t))

        print("Plotting disks")
        gg.plot_gdisks(MM,
                title="Gershgorin disks for {} - {}".format(name,method),
                filename="../img/disks/{}_{}.png".format(name,method))

        print("Performing iteration")

    # Determine error function
    # If true solution is provided, use that
    # Otherwise, check against previous iteration
    if true_sol == None:
        err_fun = lambda x: la.norm(x[:,1] - x[:,0])
    else:
        err_fun = lambda x: la.norm(x[:,1] - true_sol)

    # Initialize solution vector
    # First column is previous iteration
    # Second column is next iteration
    xx = np.zeros([nn,2])

    # Linear stationary iterative method is
    # x_{i+1} = inv(N)*P*x_i + inv(N)*b

    # Intermediate vectors:
    # C = P * x_i
    # D = inv(N) * b
    DD = tri_solve(NN,bb)

    #  Then, iterate formula is
    # x_{i+1} = inv(N)*C + D

    # Assume solution does not converge until proven otherwise
    status = 1

    # Allocate error array
    err = np.zeros(maxiter)
    niter = maxiter

    # Loop through iterations until err<tol or maxiter is reached
    for it in range(maxiter):

        # Iterate
        CC = PP @ xx[:,0]
        xx[:,1] = tri_solve(NN,CC) + DD

        # Check error
        err[it] = err_fun(xx)

        # Check convergence criterion
        if err[it] < tol:
            status = 0
            niter = it
            break

        # Update
        xx[:,0] = xx[:,1]

    return (xx, err[:niter+1], status)

# Define particular N_fun for iterative methods
N_funs = dict(
    # Jacobi
    jac  = lambda D, L, U: D ,
    # Gauss-Seidel
    gs   = lambda D, L, U: D + L ,
    # Successive over-relaxation
    sor  = lambda D, L, U, w: 1/w * D + L ,
    # Simultaneous over-relaxation
    jor  = lambda D, L, U, w: 1/w * D ,
    # Richardson
    rich = lambda D, L, U, p: -1/p * sp.eye(D.shape[0]) )

# Define callable iterative solvers
it_method = {}
for method in N_funs.keys():
    # i=i forces "early binding"
    # global variable i is saved on definition
    # Default behavior is to lookup i on function call
    # http://stackoverflow.com/questions/3431676/creating-functions-in-a-loop/
    def tmpf(AA,bb,method=method,*args,**kwargs):
        "Callable iterative solver for specified method"

        return st_itr_tmpl(AA,bb,N_funs[method],method=method,*args,**kwargs)

    # Save function to dictionary w/ appropriate name
    it_method[method] = tmpf

