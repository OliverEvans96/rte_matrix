# File Name: gen_matrix.py
# Description: Generate matrix from RTE & create image to show structure
# Created: Sun Apr 09, 2017 | 01:57pm EDT
# Last Modified: Mon Apr 24, 2017 | 10:23pm EDT

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
import numpy.polynomial.legendre as leg
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import sparse, misc, io
import scipy.sparse
import scipy.misc
import scipy.io
import IPython

# Assume space is in rescaled dimensions
# x \in [0,1), y \in [0,1)
# Kelp grows on the rope at x = 0.5 evenly in both directions
class KelpScenario(object):
    def __init__(self,surf_bc_fun,iops):
        self._grid_defined = False
        self.set_surf_bc_fun(surf_bc_fun)
        self.set_iops(*iops)

    def set_kelp(self,kelp_lengths,ind):
        self.set_ind(ind)
        self.set_kelp_lengths(kelp_lengths)
        self.set_kelp_sigma(0)

    # Actually, this won't work with angular grid defined by Legendre roots
    # (non-constant dth)
    def set_grid_spacing(self,dx,dy,dth):
        # # Set grid spacing
        # self._dx = dx
        # self._dy = dy
        # self._dth = dth

        # # Calculate number of grid points
        # self._nx = int(np.floor(1/dx))
        # self._ny = int(np.floor(1/dy))
        # self._nth = int(np.floor(2*np.pi/dth))

        # # Calculate grid points
        # self._calculate_grid()

        raise NotImplementedError

    def set_num_grid_points(self,nx,ny,nth):
        # Set number of grid points
        self._nx = nx
        self._ny = ny
        self._nth = nth

        # Calculate grid spacing
        self._dx = 1/nx
        self._dy = 1/ny
        self._dth = 2*np.pi/nth

        # Calculate grid points
        self._calculate_grid()

    # Calculate grid points after grid has been set
    def _calculate_grid(self):
        # Ensure grid has already been defined
        #assert(self._grid_defined,"Call set_grid_spacing or set_grid_divisions first")

        # Collect variable lengths for easy access
        self._var_lengths = np.array([self._nx,self._ny,self._nth],
                dtype=int)
        self._n_vars = len(self._var_lengths)

        # Define x, y and theta arrays
        # x and theta periodic, so should not contain last element
        self._xx = np.linspace(0,1,self._nx+1)[:-1]
        self._yy = np.linspace(0,1,self._ny)

        # Define theta angles according to Legendre-Gauss quadrature
        # Simulatneously define quadrature weights which will be used in scattering integral
        theta_unscaled,weights_unscaled = leg.leggauss(self._nth)

        # Scale theta from [-1,1) to [0,2pi)
        # And scale weights by pi according to change of variables
        self._theta = (theta_unscaled + 1) * np.pi
        self._lg_weights = weights_unscaled * np.pi

        # Or use evenly spaced grid
        #self._theta = np.linspace(0,2*np.pi,self._nth+1)[:-1]

        # Center of x variable
        self._xx_center = 0.5

    # Number of individuals per superindividual
    def set_ind(self,ind):
        self._ind = ind

    # Length distribution of kelp
    def set_kelp_lengths(self,kelp_lengths):
        self._kelp_lengths = kelp_lengths

    def set_kelp_sigma(self,kelp_sigma):
        self._kelp_sigma = kelp_sigma

    # Callable function of 1 variable (angle)
    def set_surf_bc_fun(self,surf_bc_fun):
        self._surf_bc_fun = surf_bc_fun

    # Set Inherent optical properties (IOPs)
    # volume scattering function (vsf)
    # absorption coefficient (abs_water/abs_kelp)
    # scattering coefficient (sct_water/abs_kelp)
    # total attenuation coefficient (atn_water/kelp)
    def set_iops(self,vsf,abs_water,sct_water,abs_kelp,sct_kelp):
        self._vsf = vsf
        self._abs_water = abs_water
        self._sct_water = sct_water
        self._abs_kelp = abs_kelp
        self._sct_kelp = sct_kelp
        self._atn_water = abs_water + sct_water
        self._atn_kelp = abs_kelp + sct_kelp

    def set_var_order(self,var_order):
        self._var_order = var_order

    # Generate 2D array of kelp density over space
    # Assume all individuals are identical. Kelp density is either ind(j)/2 or 0
    # Probability of kelp is defined as
    # P_k(i,j) = {ind(j)/2, abs(x_j - .5) <= L(j); 0, otherwise}
    def calculate_pk(self):
        self._p_k = np.zeros([self._nx,self._ny])

        # Loop over depths
        for jj in range(self._ny):
            # Loop over length
            for ii in range(self._nx):
                # If we're in the situation with no probability distribution
                if(self._kelp_sigma == 0):
                    # Nonzero only closer to center than kelp length
                    if (abs(self._xx[ii] - self._xx_center) 
                        < self._kelp_lengths[jj]):
                        self._p_k[ii,jj] = self._ind[jj] / 2
                    else:
                        self._p_k[ii,jj] = 0
                else:
                    raise NotImplementedError(
                            "Prob. dist. on lengths not supported.")

    def plot_kelp(self,imgfile=None):
        plt.clf()
        plt.imshow(self._p_k.T,extent=[0,1,0,1],
                interpolation='none')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Kelp density')
        plt.colorbar()

        if imgfile != None:
            plt.savefig(imgfile)

    def calculate_rte_matrix(self,var_order=[0,1,2]):
        # RTE matrix in linked list format
        mat = sparse.lil_matrix((*[np.prod(self._var_lengths)]*2,))
        rhs = sparse.lil_matrix((np.prod(self._var_lengths),1))
        self.set_var_order(var_order)
        
        # Alias for index function
        indx = lambda ii,jj,kk: self._matrix_index(var_order,(ii,jj,kk))

        # Loop through i, j and k in specified order
        for ind1 in range(self._var_lengths[var_order[0]]):
            # Report status
            pcnt = (ind1+1) / self._var_lengths[var_order[0]]
            print("ind1 = {}/{}: {:.1f}%".format(ind1+1,self._var_lengths[var_order[0]],100*pcnt))
            for ind2 in range(self._var_lengths[var_order[1]]):
                for ind3 in range(self._var_lengths[var_order[2]]):
                    # Extract ii, jj, and kk
                    inds = [ind1,ind2,ind3]
                    ii = inds[var_order.index(0)]
                    jj = inds[var_order.index(1)]
                    kk = inds[var_order.index(2)]
                    th = self._theta[kk]
                    row = indx(ii,jj,kk)

                    ## Cases for y derivative ##

                    # Surface BC affects downwelling radiance
                    # theta < pi (use kk < nth/2 since first endpoint
                    # is included, last is excluded
                    if jj == 0:
                        # Surface boundary condition
                        if kk < self._nth/2:
                            mat[row,row] = 1
                            rhs[row] = self._surf_bc_fun(self._theta[kk])
                            pde_flag = False
                        else:
                            # FD2 for 1st derivative in y
                            mat[row,indx(ii,jj,kk)] = -3*np.sin(th)/(2*self._dy)
                            mat[row,indx(ii,jj+1,kk)] = 2*np.sin(th)/(self._dy)
                            mat[row,indx(ii,jj+2,kk)] = -np.sin(th)/(2*self._dy)
                            pde_flag = True
                    # Bottom BC: no upwelling radiance
                    # theta >= pi
                    elif jj == self._ny-1:
                        if kk >= self._nth/2:
                            # Bottom BC
                            mat[row,row] = 1
                            rhs[row] = 0
                            pde_flag = False
                        else:
                            # BD2 for 1st derivative in y
                            mat[row,indx(ii,jj-2,kk)] = np.sin(th)/(2*self._dy)
                            mat[row,indx(ii,jj-1,kk)] = -2*np.sin(th)/(self._dy)
                            mat[row,indx(ii,jj,kk)] = 3*np.sin(th)/(2*self._dy)
                            pde_flag = True

                    # Interior of domain
                    # CD2 for y derivative
                    else:
                        # Modulus for periodic x BC
                        ip1 = np.mod(ii+1,self._nx)
                        im1 = np.mod(ii-1,self._nx)

                        # x derivative
                        mat[row,indx(ip1,jj,kk)] = np.cos(th)/(2*self._dx)
                        mat[row,indx(im1,jj,kk)] = -np.cos(th)/(2*self._dx)
                        # y derivative
                        mat[row,indx(ii,jj+1,kk)] = np.sin(th)/(2*self._dy)
                        mat[row,indx(ii,jj-1,kk)] = -np.sin(th)/(2*self._dy)

                        pde_flag = True

                    # Apply PDE when appropriate
                    # atten. & x derivative & scattering
                    if pde_flag:
                        atn_coef = (self._atn_kelp * self._p_k[ii,jj]
                                + self._atn_water * (1 - self._p_k[ii,jj]))
                        sct_coef = (self._sct_kelp * self._p_k[ii,jj]
                                + self._sct_water * (1 - self._p_k[ii,jj]))
                        # attenuation
                        mat[row,row] = atn_coef
                        # scattering
                        for ll in range(self._nth):
                            # Only consider ll != kk
                            if(ll == kk):
                                continue

                            # Theta prime - integration variable
                            thp = self._theta[ll]

                            # Take minimum distance when comparing
                            # thp and the image of thp in
                            # [2pi,4pi) and [-2pi,0)
                            # See: http://stackoverflow.com/questions/9505862/shortest-distance-between-two-degree-marks-on-a-circle
                            angle_diff = np.pi - np.abs(np.abs(
                                th - thp) - np.pi)

                            # Set scattering coefficient
                            mat[row,indx(ii,jj,ll)] = (
                                  sct_coef 
                                * self._lg_weights[ll]
                                * self._vsf(angle_diff))

        self._rte_matrix = sparse.csr_matrix(mat)
        self._rte_rhs = sparse.csr_matrix(rhs)

    # Given the order of variables in the matrix from outer to inner,
    # convert ii, jj, kk to matrix row/column #.
    # var_order is a permutation of [0,1,2] indicating the order of variables
    # 0 = xx, 1 = yy, 2 = theta
    def _matrix_index(self,var_order,indices):
        index = 0
        for nn in range(self._n_vars):
            index += (indices[var_order[nn]]
                   * np.prod(self._var_lengths[[var_order[nn+1:]]]))

        return int(index)

    def plot_rte_matrix(self,imgfile=None):
        plt.clf()
        plt.spy(self._rte_matrix,precision=0.001,markersize=.1)
        plt.title('Sparsity plot {}x{}x{}_{}{}{}'
                .format(self._nx,self._ny,self._nth,*self._var_order))

        if imgfile != None:
            plt.savefig(imgfile)

    def write_rte_matrix_png(self,imgfile):
        mat = self._rte_matrix
        mat_scaled = mat / abs(mat).max()
        mat_int = 256-abs(mat_scaled * 256).floor().astype(int).toarray()
        misc.imsave(imgfile,mat_int)

    # Like write_rte_matrix_png, but only black/white.
    # No gray to represent magnitude.
    def write_int_matrix_png(self,imgfile,guides=True):
        mat = sparse.csr_matrix(self._rte_matrix)
        mat[self._rte_matrix.nonzero()] = 256

        mat_int = abs(mat).floor().astype(int).toarray()

        if guides:
            mat_int = self._add_guides(mat_int,256)

        # Reverse colors black - white
        mat_rev = 256 - mat_int

        misc.imsave(imgfile,mat_rev)

    # Add horizontal & vertical lines of .5 & .25 of maxval at grid & subgrid
    # Grid occupies 2 pixels, subgrid occupies 1 pixel
    # Guides are added on first row & column of each block
    # (ii=1 & ii=1, jj=0 & jj=1 for grid) (ii=0,jj=0 for subgrid)
    # Only fill empty (=0) entries
    def _add_guides(self,mat,maxval=1):
        var_lengths = self._var_lengths[[self._var_order]]

        for ii in range(var_lengths[0]):
            # Grid
            ind1 = ii * np.prod(var_lengths[1:])
            # Vertical
            mat[mat[:,ind1]==0,ind1] = 0.5 * maxval
            mat[mat[:,ind1+1]==0,ind1+1] = 0.5 * maxval
            # Horizontal
            mat[ind1,mat[ind1,:]==0] = 0.5 * maxval
            mat[ind1+1,mat[ind1+1,:]==0] = 0.5 * maxval

            # Subgrid
            for jj in range(var_lengths[1]):
                ind2 = ind1 + jj * var_lengths[2]
                # Vertical
                mat[mat[:,ind2]==0,ind2] = 0.25 * maxval
                # Horizontal
                mat[ind2,mat[ind2,:]==0] = 0.25 * maxval

        return mat

    def write_rte_matrix_txt(self,out_file):
        raise NotImplementedError

    def write_rte_matrix_hdf(self,out_file):
        raise NotImplementedError

    # Allow for other variables to be stored via **kwargs
    # Union of two dicts: http://stackoverflow.com/questions/38987/how-to-merge-two-python-dictionaries-in-a-single-expression
    def write_rte_system_mat(self,out_file,**kwargs):
        io.savemat(out_file,
                {'A':self._rte_matrix,
                 'b':self._rte_rhs,
                 'x':self._rte_sol,
                 'rad':self._rad,
                 'irrad':self._irrad,
                 **kwargs}) # This trick only works for  Python3.5 +    

    def load_rte_system_mat(self,in_file,var_order=[0,1,2]):
        dct = io.loadmat(in_file)
        self._rte_matrix = dct['A']
        self._rte_rhs = dct['b']
        self.set_var_order(var_order)

    # Solve matrix equation using scipy's sparse solver.
    def solve_system(self):
        # Solve
        self._rte_sol = sparse.linalg.spsolve(self._rte_matrix,self._rte_rhs)

        # Convert to 3D array, maintaining var order
        sol3d = self._rte_sol.reshape(self._var_lengths[[self._var_order]])

        # Swap axes to order axes as [xx,yy,theta]
        self._rad = sol3d.transpose(np.argsort(self._var_order))

    # Integrate radiance over angles to get irradiance
    # Note: tensordot(*,1) - multiply & sum over last dimension of radiance (theta)
    def calc_irrad(self):
        self._irrad = np.tensordot(self._rad,self._lg_weights, axes=1)

    # Plot irradiance
    def plot_irrad(self,imgfile=None):
        plt.clf()
        plt.imshow(self._irrad.T,extent=[0,1,0,1],
                interpolation='none',cmap=cm.viridis)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Irradiance')
        plt.colorbar()

        if imgfile != None:
            plt.savefig(imgfile)

    def get_irrad(self):
        return self._irrad

