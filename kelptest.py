# File Name: kelptest.py
# Description: Generate a few matrices with gen_matrix_2d.py
# Created: Mon Apr 10, 2017 | 10:00am EDT
# Last Modified: Tue Apr 11, 2017 | 08:00pm EDT

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
import gen_matrix_2d as gm2
import itertools as it

def surf_bc_fun(th):
    return np.sin(th)

# Normalized mock VSF
def vsf(th):
    return 2 * np.exp(-th/2) / (1 - np.exp(-np.pi/2))

dx = 5e-2
dy = 5e-2
dth = np.pi/4
mesh = [dx,dy,dth]

nx = int(np.floor(1/dx))
ny = int(np.floor(1/dy))
nth = int(np.floor(2*np.pi/dth))

#kelp_lengths = np.exp(-1-np.arange(ny)/10)
kelp_lengths = np.ones(ny)
kelp_lengths[3] = 0
ind = np.ones(ny)

abs_water = 1
sct_water = 1
abs_kelp = 2
sct_kelp = 1
iops = [vsf,abs_water,sct_water,abs_kelp,sct_kelp]

scenario = gm2.KelpScenario(mesh,kelp_lengths,ind,surf_bc_fun,iops)

# What to do
gen_sparsity_plots = False
plot_irrad = True

if gen_sparsity_plots:
    # Loop through all possible variable orderings
    for ii,var_order in enumerate(it.permutations(range(3))):
        print("ii={}: {}".format(ii,var_order))
        scenario.calculate_rte_matrix(var_order)
        scenario.write_rte_matrix_png('img/rte2d_{}{}{}.png'.format(*var_order))
        scenario.write_int_matrix_png('img/int2d_{}{}{}.png'.format(*var_order))
        scenario.plot_rte_matrix('img/spy2d_{}{}{}.png'.format(*var_order))

if plot_irrad:
    print("Creating matrix")
    scenario.calculate_rte_matrix()
    #scenario.write_int_matrix_png('solve/sparsity.png')
    print("Solving system")
    scenario.solve_system()
    print("Calculating irradiance")
    scenario.calc_irrad()
    scenario.plot_irrad('solve/irrad.png')
    print("Done!")

