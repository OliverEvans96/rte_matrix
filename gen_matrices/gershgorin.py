# coding: utf-8
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
#get_ipython().magic('matplotlib inline')


# Plot gershgorin disks in the complex plane, which contain all eigenvalues of AA
# Assume AA is square
# https://en.wikipedia.org/wiki/Gershgorin_circle_theorem
# mgin is % of total width & height that each margin occupies (0<=mgin<0.5)
def plot_gdisks(AA,title=None,filename=None,fill=True,mgin=.05):
    # Dimension of matrix
    nn = AA.shape[0]
    
    # Extract diagonal elements
    DD = AA.diagonal().reshape(nn,1)
    
    # Sum of abs. in each row
    rho = np.abs(AA).sum(axis=1)
    
    # All eigenvals have mag. <= max(rho)
    max_mag = rho.max()

    print("lam <= {:.3e}".format(max_mag))
    
    # Gershorin radius = sum of abs. of non-diagonal elements
    rr = rho - np.abs(DD)

    # Set up figure
    fig = plt.figure(figsize=[10,8])
    ax = fig.gca()
    ax.set_aspect(1)
    
    # Set margin
    plt_lim = max_mag / (1 - 2*mgin)
    ax.set_xlim([-plt_lim,plt_lim])
    ax.set_ylim([-plt_lim,plt_lim])
    ax.set_xlabel(r'Re($\lambda$)')
    ax.set_ylabel(r'Im($\lambda$)')
    ax.set_title('Gershgorin Disks (n={})'.format(nn))
    
    # Plot circle containing all Gershgorin disks
    circ = plt.Circle([0,0],max_mag,fill=False,color='r')
    ax.add_artist(circ)
    
    # Plot circles in a loop (one per row)
    for ii in range(nn):
        # Center of circle in complex plane is diagonal entry
        # Radius is Gershgorin radius
        cent = [DD.real[ii,0],DD.imag[ii,0]]
        circ = plt.Circle(cent,rr[ii],
                          fill=fill,
                          color='b')
        ax.add_artist(circ)

    if filename != None:
        plt.savefig(filename)


# Plot diagonal dominance quantity (D - Q) as a function of row number
def plot_dd(AA):
    nn = AA.shape[0]
    # Diagonal dominance quantity, D - Q (pre-shift)
    ddq = np.zeros(nn)
    row_list = np.arange(nn)
    for k in row_list:
        if k%100 == 0:
            print('k={}'.format(k))
        D = np.abs(AA[k,k])
        Q = np.abs(AA[k,:k]).sum() + np.abs(AA[k,k+1:]).sum()
        ddq[k] = D - Q
    
    plt.figure(figsize=[10,6])
    plt.plot(row_list,ddq,label="$D' - Q$ (post)")
    plt.legend()
    
    print("Min (pre) : {:.3e}".format(min(ddq)))


# Given a matrix AA, return the matrix GG whose spectral radius
# determines the convergence of J, GS, and SOR methods.
# A = N - P
# G = N^{-1}P
# w is SOR parameter
def conv_mat(AA,w=1.5):
    # Determine whether to use standard or sparse functions
    if sp.issparse(AA):
        print("Sparse")
        # scipy.sparse
        pkg = sp
        # Create diagonal matrix from vector of entries
        diag = sp.diags
        # Which algorithm to use to solve a triangular system
        tri_solve = sla.spsolve
    else:
        print("Dense")
        # numpy
        pkg = np
        diag = np.diag
        tri_solve = lp.trtrs
        
    # NN, PP, and GG are lists of arrays (one for each method)
    # 0 = Jacobi,
    # 1 = Gauss-Seidel,
    # 2 = SOR
    
    # Number of iterative methods to analyze
    nim = 3

    # Allocate blank lists
    NN = [0] * nim
    PP = [0] * nim
    GG = [0] * nim
    
    # Isolate diagonal & strictly lower & upper triangular parts
    DD = diag(AA.diagonal())
    LL = pkg.tril(AA) - DD
    UU = pkg.triu(AA) - DD
    
    # lp.trtrs solves a triangular system
    # http://www.netlib.org/lapack/explore-html/da/dba/group__double_o_t_h_e_rcomputational_ga4e87e579d3e1a56b405d572f868cd9a1.html#ga4e87e579d3e1a56b405d572f868cd9a1
    
    # Jacobi
    NN[0] = DD
    #PP[0] = -(LL + UU)
    
    # Gauss-Seidel
    NN[1] = DD + LL
    #PP[1] = -UU
    
    # SOR
    NN[2] = 1/w * DD + LL
    #PP[2] = (1/w - 1) * DD - UU
    
    # Calculate GG
    for ii in range(nim):
        print("i={}".format(ii))
        # Calculate PP
        PP[ii] = NN[ii] - AA

        # Convert to CSC format for efficient solution
        if sp.issparse(AA):
            NN[ii] = NN[ii].tocsc()
            PP[ii] = PP[ii].tocsc()

        GG[ii] = tri_solve(NN[ii],PP[ii])
    
    return GG    


# Load matrix
#print("Load matrix")
#A1,b1 = [io.loadmat('../mat/ddom_20x20x24_012.mat')[key] for key in ['A','b']]
#
#print("Running conv_mat")
#GG = conv_mat(A1)
#
#print("Plot disks")
#names = ['jac','gs','sor']
#for ii in range(3):
#    print("ii={}".format(ii))
#    plot_gdisks(GG[ii],
#            title=names[ii],
#            filename=names[ii]+'.png',
#            fill=False)
#
#print("Done")
#IPython.embed()
