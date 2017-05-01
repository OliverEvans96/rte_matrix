# File Name: LU_vs_LUspares.py
# Description: plot computation time vs matrix size for LU.py and LUsparse.py
# Author: Christopher Parker
# Created: Sun Apr 30, 2017 | 07:10P EDT
# Last Modified: Sun Apr 30, 2017 | 07:15P EDT

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
from LU import LU
from LUsparse import LUsparse

LU_answers = [LU('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_10x10x16_012.mat'), LU('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_20x20x12_012.mat'), LU('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_20x20x24_012.mat'), LU('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_20x20x32_012.mat')]#, LU('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_50x50x16_012.mat'), LU('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_50x50x24_012.mat'), LU('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_50x50x32_012.mat')]

LU_times = [LU_answers[0][1], LU_answers[1][1], LU_answers[2][1], LU_answers[3][1]]#, LU_answers[4][1], LU_answers[5][1], LU_answers[6][1]]


LUsparse_answers = [LUsparse('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_10x10x16_012.mat'), LUsparse('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_20x20x12_012.mat'), LUsparse('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_20x20x24_012.mat'), LUsparse('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_20x20x32_012.mat')]#, LUsparse('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_50x50x16_012.mat'), LUsparse('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_50x50x24_012.mat'), LUsparse('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_50x50x32_012.mat')]

LUsparse_times = [LUsparse_answers[0][1], LUsparse_answers[1][1], LUsparse_answers[2][1], LUsparse_answers[3][1]]#, LUsparse_answers[4][1]], LUsparse_answers[5][1], LUsparse_answers[6][1]]

matrix_sizes = [len(LUsparse_answers[0][0]), len(LUsparse_answers[1][0]), len(LUsparse_answers[2][0]), len(LUsparse_answers[3][0])]#, len(LUsparse_answers[4][0]), len(LUsparse_answers[5][0]), len(LUsparse_answers[6][0])]

plt.plot(matrix_sizes,LU_times,'k')
plt.plot(matrix_sizes,LUsparse_times,'b')

plt.xlabel('Matrix Size')
plt.ylabel('Computation Time (s)')

plt.savefig('LU_vs_LUsparse.png')
plt.show()
