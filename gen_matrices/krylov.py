# File Name: krylov.py
# Description: Krylov subspace (non-stationalry iterative) methods
# Created: Thu May 04, 2017 | 01:18am EDT
# Last Modified: Thu May 04, 2017 | 07:17am EDT

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

import scipy.sparse.linalg as sla
import time

# Time function & return beginning & ending time
def time_func(f):
    def time_wrapper(*args,**kwargs):
        t1 = time.time()
        result = f(*args,**kwargs)
        t2 = time.time()

        # Return expanded arg list if appropriate
        if type(result) in (tuple,list):
            return (*result,t1,t2)
        else:
            return (result,t1,t2)

    return time_wrapper

krylov_method = dict(
        bicg     = time_func(sla.bicg),
        bicgstab = time_func(sla.bicgstab),
        cgs      = time_func(sla.cgs),
        gmres    = time_func(sla.gmres),
        lgmres   = time_func(sla.lgmres),
        qmr      = time_func(sla.qmr),
        lsqr     = time_func(sla.lsqr),
        lsmr     = time_func(sla.lsmr))

