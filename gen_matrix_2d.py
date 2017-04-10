# File Name: gen_matrix.py
# Description: Generate matrix from RTE & create image to show structure
# Created: Sun Apr 09, 2017 | 01:57pm EDT
# Last Modified: Mon Apr 10, 2017 | 03:06pm EDT

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
import matplotlib.pyplot as plt
from scipy import sparse, misc
import IPython

# Assume space is in rescaled dimensions
# x \in [0,1), y \in [0,1)
# Kelp grows on the rope at x = 0.5 evenly in both directions
class KelpScenario(object):
    def __init__(self,mesh,kelp_lengths,ind,surf_bc_fun,iops):
        self.set_spatial_mesh(*mesh)
        self.set_ind(ind)
        self.set_kelp_lengths(kelp_lengths)
        self.set_kelp_sigma(0)
        self.set_surf_bc_fun(surf_bc_fun)
        self.set_iops(*iops)

        self.calculate_pk()

    def set_spatial_mesh(self,dx,dy,dth):
        # Set grid size
        self._dx = dx
        self._dy = dy
        self._dth = dth

        # Calculate number of grid points
        self._nx = int(np.floor(1/dx))
        self._ny = int(np.floor(1/dy))
        self._nth = int(np.floor(2*np.pi/dth))
        self._var_lengths = np.array([self._nx,self._ny,self._nth],
                dtype=int)
        self._n_vars = len(self._var_lengths)

        # Define x, y and theta arrays
        # x and theta periodic, so should not contain last element
        self._xx = np.linspace(0,1,self._nx+1)[:-1]
        self._yy = np.linspace(0,1,self._ny)
        self._theta = np.linspace(0,2*np.pi,self._nth+1)[:-1]

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
    # absorption coefficient (abs_coef)
    # scattering coefficient (sct_coef)
    # total attenuation coefficient (atn_coef)
    def set_iops(self,vsf,abs_coef,sct_coef):
        self._vsf = vsf
        self._abs_coef = abs_coef
        self._sct_coef = sct_coef
        self._atn_coef = abs_coef + sct_coef

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
                        <= self._kelp_lengths[jj]):
                        self._p_k[ii,jj] = self._ind[jj] / 2
                    else:
                        self._p_k[ii,jj] = 0
                else:
                    raise NotImplementedError(
                            "Prob. dist. on lengths not supported.")

    def calculate_rte_matrix(self,var_order):
        # RTE matrix in linked list format
        mat = sparse.lil_matrix((*[np.prod(self._var_lengths)]*2,))
        rhs = sparse.lil_matrix((np.prod(self._var_lengths),1))
        
        # Alias for index function
        indx = lambda ii,jj,kk: self._matrix_index(var_order,(ii,jj,kk))

        # Loop through i, j and k in specified order
        for ind1 in range(self._var_lengths[var_order[0]]):
            for ind2 in range(self._var_lengths[var_order[1]]):
                for ind3 in range(self._var_lengths[var_order[2]]):
                    # Extract ii, jj, and kk
                    inds = [ind1,ind2,ind3]
                    ii = inds[var_order.index(0)]
                    jj = inds[var_order.index(1)]
                    kk = inds[var_order.index(2)]
                    th = self._theta[kk]
                    row = indx(ii,jj,kk)
                    #print("(ii,jj,kk) = ({},{},{})".format(ii,jj,kk))

                    # Surface BC affects downwelling radiance
                    # theta < pi (use kk <= nth/2 since theta[nth/2] < pi)
                    if jj == 0:
                        if kk <= self._nth/2:
                            mat[row,row] = 1
                            rhs[row] = self._surf_bc_fun(self._theta[kk])
                    # Bottom BC: no upwelling radiance
                    # theta >= pi
                    elif jj == self._ny-1:
                        if kk > self._nth/2:
                            mat[row,row] = 1
                            rhs[row] = 0
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

                        # attenuation
                        mat[row,row] = self._atn_coef
                        # scattering
                        for ll in range(self._nth):
                            # Only consider ll != kk
                            if(ll == kk):
                                continue

                            # Theta prime - integration variable
                            thp = self._theta[ll]

                            # Shortest distance in periodic var.
                            # See: http://stackoverflow.com/questions/9505862/shortest-distance-between-two-degree-marks-on-a-circle
                            angle_diff = np.pi - np.abs(np.abs(
                                th - thp) - np.pi)

                            # Take minimum distance when comparing
                            # thp and the image of thp in
                            # [2pi,4pi) and [-2pi,0)
                            mat[row,indx(ii,jj,ll)] = (
                                  self._sct_coef 
                                * self._dth
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
        if imgfile != None:
            plt.savefig(imgfile)

    def write_rte_matrix_png(self,imgfile):
        mat = self._rte_matrix
        mat_scaled = mat / abs(mat).max()
        mat_int = 256-abs(mat_scaled * 256).floor().astype(int).toarray()
        misc.imsave(imgfile,mat_int)

    # Like write_rte_matrix_png, but only black/white.
    # No gray to represent magnitude.
    def write_int_matrix_png(self,imgfile):
        mat = sparse.csr_matrix(self._rte_matrix)
        mat[self._rte_matrix.nonzero()] = 256
        mat_int = 256-abs(mat).floor().astype(int).toarray()
        misc.imsave(imgfile,mat_int)

    def write_rte_matrix_txt(self,out_file):
        pass

    def write_rte_matrix_hdf(self,out_file):
        pass
