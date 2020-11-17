# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 15:27:37 2019

@author: yadav
"""

#code to obtain normalization constant of density

import math
import cmath
import numpy
#import contour_integral
import integration_trapezoidal_nonuniformspacing

gamma=0.4
data1 = numpy.loadtxt("input_output_files/rho_2.0/gamma_"+str(gamma)+"/density/renormalized_density_psi_method_2_epsi=1e-4_gamma="+str(gamma)+"_theta=2.0_rho=2.0_18000points_corrected4_iter23.0.txt",float)

f_out2=file("input_output_files/rho_2.0/gamma_"+str(gamma)+"/density/renormalized_density_psi_method_2_epsi=1e-4_gamma="+str(gamma)+"_theta=2.0_rho=2.0_18000points_corrected4_iter23.0_normalized.txt","w")
f_out=file("input_output_files/rho_2.0/density/normalization_constant_density_corrected4_iter10.0.txt","a")

#following line adds header information to output data file. 
#f_out.write("this file has normalization constant for density corresponding to given gamma. 1st column-gamma, 2nd column-normalization constant."+'\n')    # comment this line after one-time use

x = data1[:, 0]    # Array Slicing: Accessing array's first column(i.e. index-zero column), https://jakevdp.github.io/PythonDataScienceHandbook/02.02-the-basics-of-numpy-arrays.html
f_x = data1[:, 1]    # Array Slicing: Accessing array's first column(i.e. index-zero column)

N0 = abs(integration_trapezoidal_nonuniformspacing.int_trap_non_uniform(x,f_x))
print N0

f_out.write(str(gamma)+" "+str(N0)+'\n')
f_out.close()    # () at the end is necessary to close the file

for i in range(len(data1)):    # function len() on array gives no. of rows of array 
    x1 = data1[i,0]
    normal_density = (1.0/N0)*data1[i,1]
    f_out2.write(str(x1)+" "+str(normal_density)+'\n')    
#    print int(i)
    
f_out2.close()    # () at the end is necessary to close the file         