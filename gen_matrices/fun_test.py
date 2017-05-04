# File Name: fun_test.py
# Description: Check passing functions
# Author: Oliver Evans
# Created: Wed May 03, 2017 | 02:19pm EDT
# Last Modified: Wed May 03, 2017 | 03:29pm EDT

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

import copy

n=3

def mult(a,b):
    return a*b

print("Define")
fun_list = {}
for i in range(n):
    print("f[{}](x) = {}*x".format(i,i))
    # i=i forces "early binding"
    # i is saved on definition, not on function call
    # http://stackoverflow.com/questions/3431676/creating-functions-in-a-loop/
    def tmpf(x,i=i):
        return mult(x,i)
    fun_list[i] = tmpf

print()
print("Call")
x=2
print("x={}".format(x))
for i in range(n):
    print("f[{}]({}) = {}".format(i,x,fun_list[i](x)))

i=0
print("f[{}]({}) = {}".format(i,x,fun_list[i](x)))
