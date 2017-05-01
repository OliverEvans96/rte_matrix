# File Name: QR_vs_SVD_vs_LU.py
# Description: 
# Author: Christopher Parker
# Created: Sun Apr 30, 2017 | 10:32P EDT
# Last Modified: Sun Apr 30, 2017 | 11:58P EDT

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
from QR import QR
from SVD import SVD

LU_answers = [LU('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_10x10x16_012.mat'), LU('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_20x20x12_012.mat'), LU('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_20x20x24_012.mat'), LU('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_20x20x32_012.mat')]#, LU('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_50x50x16_012.mat'), LU('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_50x50x24_012.mat'), LU('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_50x50x32_012.mat')]

LU_times = [LU_answers[0][1], LU_answers[1][1], LU_answers[2][1], LU_answers[3][1]]#, LU_answers[4][1], LU_answers[5][1], LU_answers[6][1]]


QR_answers = [QR('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_10x10x16_012.mat'), QR('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_20x20x12_012.mat'), QR('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_20x20x24_012.mat'), QR('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_20x20x32_012.mat')]#, QR('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_50x50x16_012.mat'), QR('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_50x50x24_012.mat'), QR('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_50x50x32_012.mat')]

QR_times = [QR_answers[0][1], QR_answers[1][1], QR_answers[2][1], QR_answers[3][1]]#, QR_answers[4][1]], QR_answers[5][1], QR_answers[6][1]]


SVD_answers = [SVD('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_10x10x16_012.mat'), SVD('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_20x20x12_012.mat'), SVD('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_20x20x24_012.mat'), SVD('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_20x20x32_012.mat')]#, SVD('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_50x50x16_012.mat'), SVD('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_50x50x24_012.mat'), SVD('/Users/christopher/Documents/anaII/rte_matrix/mat/kelp1_50x50x32_012.mat')]

SVD_times = [SVD_answers[0][1], SVD_answers[1][1], SVD_answers[2][1], SVD_answers[3][1]]#, SVD_answers[4][1]], SVD_answers[5][1], SVD_answers[6][1]]

matrix_sizes = [len(QR_answers[0][0]), len(QR_answers[1][0]), len(QR_answers[2][0]), len(QR_answers[3][0])]#, len(QR_answers[4][0]), len(QR_answers[5][0]), len(QR_answers[6][0])]

plt.plot(matrix_sizes,LU_times,'k')
plt.plot(matrix_sizes,QR_times,'b')
plt.plot(matrix_sizes,SVD_times,'r')

plt.xlabel('Matrix Size')
plt.ylabel('Computation Time (s)')

plt.savefig('QR_vs_SVD_vs_LU.png')
plt.show()
