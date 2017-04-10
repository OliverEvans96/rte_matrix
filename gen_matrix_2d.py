# File Name: gen_matrix.py
# Description: Generate matrix from RTE & create image to show structure
# Created: Sun Apr 09, 2017 | 01:57pm EDT
# Last Modified: Sun Apr 09, 2017 | 08:11pm EDT

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

from numpy import *
from matplotlib.pyplot import *

# Assume space is in rescaled dimensions
# x \in [0,1), y \in [0,1)
# Kelp grows on the rope at x = 0.5 evenly in both directions
class KelpScenario(object):
    def __init__(self):
        pass

    def set_spatial_mesh(self,dx,dy,dth):
        # Set grid size
        self._dx = dx
        self._dy = dy
        self._dth = dth

        # Calculate number of grid points
        self._nx = 1/dx
        self._ny = 1/dy
        self._nth = 1/dth

        # Define x, y and theta arrays
        self._xx = linspace(0,1,self._nx)
        self._yy = linspace(0,1,self._ny)
        self._theta = linspace(0,1,self._nth)

        # Center of x variable
        self._xx_center = 0.5

    # Number of individuals per superindividual
    def set_ind(self,ind):
        self._ind = ind

    # Length distribution of kelp
    def set_kelp_lengths(self,kelp_lengths):
        self._kelp_lengths = kelp_lengths

    def set_kelp_sigma(srlf,kelp_sigma):
        self._kelp_sigma = kelp_sigma

    # Callable function of 1 variable (angle)
    def set_surf_bc_fun(self,surf_bc_fun):
        self._surf_bc_fun = surf_bc_fun

    # Set Inherent optical properties (IOPs)
    # volume scattering function (vsf) - normalized
    # absorption coefficient (abs_coef)
    # scattering coefficient (sct_coef)
    def set_iops(self,vsf,abs_coef,sct_coef):
        self._vsf = vsf
        self._abs_coef = abs_coef
        self._sct_coef = sct_coef

    # Generate 2D array of kelp density over space
    # Assume all individuals are identical. Kelp density is either ind(j)/2 or 0
    # Probability of kelp is defined as
    # P_k(i,j) = {ind(j)/2, abs(x_j - .5) <= L(j); 0, otherwise}
    def calculate_pk(self):
        # Loop over depths
        for jj in range(self._ny):
            # Loop over length
            for ii in range(self._nx):
                # If we're in the situation with no probability distribution
                if(self._kelp_sigma == 0):
                    # Nonzero only closer to center than kelp length
                    if abs(self._xx[ii] - self._xx_center) <= kelp_lengths(jj):
                        p_k[ii,jj] = self._ind(jj) / 2
                    else:
                        p_k[ii,jj] = 0
                else:
                    raise NotImplementedError("Prob. dist. on lengths not supported.")

    def calculate_rte_matrix(self):
        pass

    def plot_rte_matrix(self,imgfile=None):
        pass

    def write_rte_matrix_txt(self,out_file):
        pass

    def write_rte_matrix_hdf(self,out_file):
        pass
