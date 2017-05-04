# File Name: run_trials.py
# Description: Try combinations of matrix size and method,
# collecting performance data for each.
# Created: Thu May 04, 2017 | 05:17am EDT
# Last Modified: Thu May 04, 2017 | 05:56am EDT

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
from scipy import io
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as sla
import glob

import krylov
import direct
import stat_it

# dictionary of dictionaries of solvers
solvers = dict(
        krylov = krylov.krylov_method, 
        direct = direct.direct_method, 
        stat_it = stat_it.stat_method)

# List of matrices
matrices = glob.glob('../mat/ddom_*_210.mat')

# Extract matrix sizes
names = [mat.split('_')[1] for mat in matrices]
dims = np.array([[int(dim) for dim in name.split('x')])
sizes = np.prod(dims, axis=1)

# Loop through 
for solver 

