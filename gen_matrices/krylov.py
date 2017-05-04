# File Name: krylov.py
# Description: Krylov subspace (non-stationalry iterative) methods
# Created: Thu May 04, 2017 | 01:18am EDT
# Last Modified: Thu May 04, 2017 | 01:38am EDT

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

krylov_method = dict(
        bicg = sla.bicg,
        bicgstab = sla.bicgstab,
        cgs = sla.cgs,
        gmres = sla.gmres,
        lgmres = sla.lgmres,
        lmres = sla.lmres,
        qmr = sla.qmr,
        lsqr = sla.lsqr,
        lsmr = sla.lsmr)

