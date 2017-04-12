# File Name: kelptest.py
# Description: Generate a few matrices with gen_matrix_2d.py
# Created: Mon Apr 10, 2017 | 10:00am EDT
# Last Modified: Wed Apr 12, 2017 | 07:54am EDT

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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import gen_matrix_2d as gm2
import itertools as it

def surf_bc_fun(th):
    return np.sin(th)

# Normalized mock VSF
def vsf(th):
    return 2 * np.exp(-th/2) / (1 - np.exp(-np.pi/2))

dx = 1e-2
dy = 1e-2
dth = np.pi/12
mesh = [dx,dy,dth]

nx = int(np.floor(1/dx))
ny = int(np.floor(1/dy))
nth = int(np.floor(2*np.pi/dth))

xx = np.linspace(0,1,nx)
yy = np.linspace(0,1,ny)
theta = np.linspace(0,2*np.pi,nth)

# Made up shape. Little kelp on top, bulge near middle, zero at bottom.
kelp_lengths = 5 * (1 - yy) ** 2 * np.exp(5 * yy - 4)
# Number of individual kelps in each depth layer - more towards the top
ind = 2 - yy

# Assume that kelp and water scatter the same,
# but kelp absorbs much more light
abs_water = 1
sct_water = 1
abs_kelp = 20
sct_kelp = 1
iops = [vsf,abs_water,sct_water,abs_kelp,sct_kelp]

scenario = gm2.KelpScenario(mesh,kelp_lengths,ind,surf_bc_fun,iops)

# What to do
gen_sparsity_plots = True
plot_kelp = False
plot_irrad = False

print("{}x{}x{}".format(nx,ny,nth))

if gen_sparsity_plots:
    # Loop through all possible variable orderings
    for ii,var_order in enumerate(it.permutations(range(3))):
        print()
        print("ii={}: {}".format(ii,var_order))

        # Determine common name for files
        # kelp1_[variable dimensions]_[variable order]
        name = ('kelp1_{}x{}x{}_{}{}{}'
                .format(nx,ny,nth,*var_order))

        print("Creating matrix")
        scenario.calculate_rte_matrix(var_order)

        print("Saving files")
        # Save mat file
        scenario.write_rte_system_mat('mat/'+name)
        # Save sparsity plots - one coarse (spy) & one precise (int)
        #scenario.write_int_matrix_png('img/sparsity/int_'+name+'.png')
        scenario.plot_rte_matrix('img/sparsity/spy_'+name+'.png')

        # Solve system & plot result
        print("Solving system")
        scenario.solve_system()
        print("Calculating irradiance")
        scenario.calc_irrad()
        scenario.plot_irrad('img/irrad/irrad_'+name+'.png')

if plot_kelp:
    print("Plotting kelp")
    plt.figure(1)
    scenario.plot_kelp('solve/kelp.png')

if plot_irrad:
    print("Creating matrix")
    plt.figure(2)
    scenario.calculate_rte_matrix()
    #scenario.write_int_matrix_png('solve/sparsity.png')
    print("Solving system")
    scenario.solve_system()
    print("Calculating irradiance")
    scenario.calc_irrad()
    scenario.plot_irrad('solve/irrad.png')
    print("Done!")

