# File Name: conv_study.py
# Description: Generate a set of matrices for RTE2D to study numerical
# convergence
# Created: Mon Apr 10, 2017 | 10:00am EDT
# Last Modified: Wed Apr 12, 2017 | 04:20pm EDT

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

import matplotlib
#matplotlib.use('Qt5Agg') # interactive
#matplotlib.use('Agg') # noninteractive
import numpy as np
import time
from scipy import io
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import gen_matrix_2d as gm2
import itertools as it
import IPython

## SCENARIO PARAMETERS ##

def surf_bc_fun(th):
    return np.sin(th)

# Normalized mock VSF
def vsf(th):
    return 2 * np.exp(-th/2) / (1 - np.exp(-np.pi/2))

# Assume that kelp and water scatter the same,
# but kelp absorbs much more light
abs_water = 1
sct_water = 1
abs_kelp = 5
sct_kelp = 1
iops = [vsf,abs_water,sct_water,abs_kelp,sct_kelp]

## CONVERGENCE STUDY ##
# Here, we're specifying number of points, including the left endpoint, for both spatial and angular grids
# Also, dx = dy.

# Factor by which to increase number of grid points in each grid refinement.
# Could be different for spatial vs. angular
# Should be an integer for simplicity
rs_ref = 2
ra_ref = 2

# Number of spatial grid trials
nsg = 3
# Number of angular grid trials
nag = 3

# Minimum & Maximum grid size
# Spatial
ns_min = 10
ns_max = ns_min * rs_ref ** nsg
# Angular
na_min = 10
na_max = na_min * ra_ref ** nag

# Array of grid sizes to loop over
# First axis is spatial, second is angular
spatial_grids = ns_min * rs_ref ** np.arange(nsg)
angular_grids = na_min * ra_ref ** np.arange(nag)

# Total error array (2D)
# 1st dimension - spatial trials
# 2nd dimension - angular trials
tot_err = np.zeros([nsg,nag])

# We can't compare radiance values among angular grid sizes since we define theta according to
# the roots of the Legendre polynomials, but we can compare irradiances.

# Compare at points present in coarsest resolution
# Compare to "true solution", i.e., highest resolution,
# So we have to calculate highest resolution grid first.

# "True irradiance"
# values from highest resolution trial
# grid points from lowest resolution trial
irrad_true = np.zeros([ns_min,ns_min])

# Create kelp scenario
scenario = gm2.KelpScenario(surf_bc_fun,iops)

print("Convergence study")

# Loop over spatial/angular grid sizes (highest first)
for ii in range(nsg-1,-1,-1):
    for jj in range(nag-1,-1,-1):

        ns = spatial_grids[ii]
        na = angular_grids[jj]

        print()
        print("ns={} x na={}".format(ns,na))
        
        # Perform calculations for this trial
        scenario.set_num_grid_points(ns,ns,na)

        ## KELP DISTRIBUTION ##
        yy = scenario._yy
        # Made up shape. Little kelp on top, bulge near middle, zero at bottom.
        kelp_lengths = 5 * (1 - yy) ** 2 * np.exp(5 * yy - 4)
        # Number of individual kelps in each depth layer - more towards the top
        ind = 2 - yy
        
        # Input the distribution to the scenario
        print("Set kelp")
        scenario.set_kelp(kelp_lengths,ind)
        # Calculate probability of kelp over 2D space
        print("Calculate pk")
        scenario.calculate_pk()

        # Generate system & solve 
        print("Gen matrix")
        scenario.calculate_rte_matrix()
    
        # Solve system - save the amount of time required to solve
        print("Solve system")
        t1 = time.time()
        scenario.solve_system()
        solve_time = time.time() - t1
        print("dt = {:.2e}".format(solve_time))

        print("Calc irrad")
        scenario.calc_irrad()
        irrad = scenario.get_irrad()

        # Length of skips in order to only access LRGPs
        skip_len = rs_ref ** ii

        # Save irradiance values from highest resolution spatio-angular grid
        # Only need to save at lowest resolution grid points
        if ( (ii == nsg-1) and (jj == nag-1) ):
            # Save it, man!
            irrad_true = irrad[::skip_len,::skip_len]
            
            # Error is zero! (since this is the "true" trial)
            tot_err[ii,jj] = 0

        # For others, we compare irradiance w/ "true" to calculate error
        else:
            # Divide by number of points = ns_min **2 so that error is average per point
            tot_err[ii,jj] = np.sum(np.abs(irrad[::skip_len,::skip_len] - irrad_true)) / ns_min ** 2

        # Report error
        print("err = {:.2e}".format(tot_err[ii,jj]))

        # Save system
        name = 'conv_{}x{}'.format(ns,na)
        out_file = '../mat/conv/' + name + '.mat'
        # Save error as well (not really necessary)
        scenario.write_rte_system_mat(out_file,
                err=tot_err[ii,jj],
                py_solve_time=solve_time)

        # Save irradiance plot
        scenario.plot_irrad('../img/irrad/conv/' + name + '.png')

# After all trials have been calculated, save tot_err
io.savemat('../mat/conv/tot_error.mat',{'tot_err': tot_err})

print("Done!")

