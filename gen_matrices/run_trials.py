# File Name: run_trials.py
# Description: Try combinations of matrix size and method,
# collecting performance data for each.
# Created: Thu May 04, 2017 | 05:17am EDT
# Last Modified: Thu May 04, 2017 | 07:31am EDT

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

import IPython
import numpy as np
from scipy import io
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as sla
import pandas as pd
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
dims = np.array([[int(dim) for dim in name.split('x')] for name in names])
sizes = np.prod(dims, axis=1)

# Loop through matrices for krylov methods
slv_class='krylov'
nmatrices = len(matrices)
nmethods = len(solvers[slv_class].keys())
start_times = np.zeros([nmatrices,nmethods])
end_times = np.zeros([len(matrices),len(solvers[slv_class].keys())])

# Pandas DataFrame for each class of methods
# Row for each trial
# 3 columns: start_time, end_time, info (failed to converge?)
data = {}
data[slv_class] = pd.DataFrame(
    np.zeros([nmatrices*nmethods,7]),
    columns=['method','mat_name','mat_size','info','start_time','end_time','dt'])

### These solvers have a callback method if we want to calculate errors ###

slv_keys = list(solvers[slv_class].keys())

for ii,matrix in enumerate(matrices):
    mat_dct = io.loadmat(matrix)
    A = mat_dct['A']
    b = mat_dct['b'].toarray()
    for jj, method in enumerate(slv_keys):
        print("i={}, j={}".format(ii,jj))
        kk = ii * nmethods + jj
        out = solvers[slv_class][method](A,b)

        data[slv_class].loc[kk,'method'] = slv_keys[jj]
        data[slv_class].loc[kk,'mat_name'] = names[ii]
        data[slv_class].loc[kk,'mat_size'] = sizes[ii]
        data[slv_class].loc[kk,'info'] = out[1]
        data[slv_class].loc[kk,'start_time'] = out[2]
        data[slv_class].loc[kk,'end_time'] = out[3]
        data[slv_class].loc[kk,'dt'] = out[3] - out[2]

io.savemat('../data/{}'.format(slv_class),data)
