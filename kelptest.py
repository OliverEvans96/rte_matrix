# File Name: kelptest.py
# Description: Generate a few matrices with gen_matrix_2d.py
# Created: Mon Apr 10, 2017 | 10:00am EDT
# Last Modified: Mon Apr 10, 2017 | 03:08pm EDT

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

def vsf(th):
    return np.exp(-th/2)

dx = 1e-1
dy = 1e-1
dth = np.pi/8
mesh = [dx,dy,dth]

nx = int(np.floor(1/dx))
ny = int(np.floor(1/dy))
nth = int(np.floor(2*np.pi/dth))

kelp_lengths = np.ones(ny)
ind = np.ones(ny)

abs_coef = 1
sct_coef = 1
iops = [vsf,abs_coef,sct_coef]

scenario = gm2.KelpScenario(mesh,kelp_lengths,ind,surf_bc_fun,iops)

# Loop through all possible variable orderings
for ii,var_order in enumerate(it.permutations(range(3))):
    print("ii={}: {}".format(ii,var_order))
    scenario.calculate_rte_matrix(var_order)
    scenario.write_rte_matrix_png('img/rte2d_{}{}{}.png'.format(*var_order))
    scenario.write_int_matrix_png('img/int2d_{}{}{}.png'.format(*var_order))
    scenario.plot_rte_matrix('img/spy2d_{}{}{}.png'.format(*var_order))
