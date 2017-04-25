# File Name: kelptest.py
# Description: Generate a few matrices with gen_matrix_2d.py
# Created: Mon Apr 10, 2017 | 10:00am EDT
# Last Modified: Mon Apr 24, 2017 | 10:18pm EDT

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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import gen_matrix_2d as gm2
import itertools as it
import IPython

def surf_bc_fun(th):
    return np.sin(th)

# Normalized mock VSF
def vsf(th):
    return 2 * np.exp(-th/2) / (1 - np.exp(-np.pi/2))

nx = 10
ny = 10
nth = 16

xx = np.linspace(0,1,nx)
yy = np.linspace(0,1,ny)

# Made up shape. Little kelp on top, bulge near middle, zero at bottom.
kelp_lengths = 5 * (1 - yy) ** 2 * np.exp(5 * yy - 4)
# Number of individual kelps in each depth layer - more towards the top
ind = 2 - yy

# Assume that kelp and water scatter the same,
# but kelp absorbs much more light
abs_water = 1
sct_water = 1
abs_kelp = 5
sct_kelp = 1
iops = [vsf,abs_water,sct_water,abs_kelp,sct_kelp]

scenario = gm2.KelpScenario(surf_bc_fun,iops)
scenario.set_kelp(kelp_lengths,ind)
scenario.set_num_grid_points(nx,ny,nth)
scenario.calculate_pk()

# What to do
gen_sparsity_plots = True
interactive_load_mat = False
plot_kelp = False
plot_irrad = True

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
        # Solve system & plot result
        print("Solving system")
        scenario.solve_system()
        print("Calculating irradiance")

        scenario.calc_irrad()
        scenario.plot_irrad('../img/irrad/irrad_'+name+'.png')

        # Save mat file
        scenario.write_rte_system_mat('../mat/'+name)
        # Save sparsity plots - one coarse (spy) & one precise (int)
        scenario.write_int_matrix_png('../img/sparsity/int_'+name+'.png')
        scenario.plot_rte_matrix('../img/sparsity/spy_'+name+'.png')

if plot_kelp:
    print("Plotting kelp")
    plt.figure(1)
    scenario.plot_kelp('../img/solve/kelp.png')

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

if interactive_load_mat:
    print("Loading matrix")
    scenario.load_rte_system_mat('mat/kelp1_50x50x32_012.mat',[0,1,2])
    print("Finished loading")
    IPython.embed()
