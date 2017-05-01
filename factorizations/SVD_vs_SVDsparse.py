# File Name: SVD_vs_SVDsparse.py
# Description: 
# Author: Christopher Parker
# Created: Sun Apr 30, 2017 | 07:10P EDT
# Last Modified: Sun Apr 30, 2017 | 08:46P EDT

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#                           GNU GPL LICENSE                            #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#                                                                      #
# Copyright Christopher Parker 2017 <cjp65@case.edu>                   #
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

import matplotlib.pyplot as plt
from SVD import SVD
from SVDsparse import SVDsparse

SVD_answers = [SVD('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_10x10x16_012.mat'), SVD('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_20x20x12_012.mat'), SVD('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_20x20x24_012.mat'), SVD('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_20x20x32_012.mat')]#, SVD('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_50x50x16_012.mat'), SVD('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_50x50x24_012.mat'), SVD('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_50x50x32_012.mat')]

SVD_times = [SVD_answers[0][1], SVD_answers[1][1], SVD_answers[2][1], SVD_answers[3][1]]#, SVD_answers[4][1], SVD_answers[5][1], SVD_answers[6][1]]


SVDsparse_answers = [SVDsparse('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_10x10x16_012.mat'), SVDsparse('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_20x20x12_012.mat'), SVDsparse('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_20x20x24_012.mat'), SVDsparse('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_20x20x32_012.mat')]#, SVDsparse('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_50x50x16_012.mat'), SVDsparse('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_50x50x24_012.mat'), SVDsparse('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_50x50x32_012.mat')]

SVDsparse_times = [SVDsparse_answers[0][1], SVDsparse_answers[1][1], SVDsparse_answers[2][1], SVDsparse_answers[3][1]]#, SVDsparse_answers[4][1]], SVDsparse_answers[5][1], SVDsparse_answers[6][1]]

matrix_sizes = [len(SVDsparse_answers[0][0]), len(SVDsparse_answers[1][0]), len(SVDsparse_answers[2][0]), len(SVDsparse_answers[3][0])]#, len(SVDsparse_answers[4][0]), len(SVDsparse_answers[5][0]), len(SVDsparse_answers[6][0])]

plt.semilogy(matrix_sizes,SVD_times,'k')
plt.semilogy(matrix_sizes,SVDsparse_times,'b')

plt.xlabel('Matrix Size')
plt.ylabel('Computation Time (s)')

plt.savefig('SVD_vs_SVDsparse.png')
plt.show()
