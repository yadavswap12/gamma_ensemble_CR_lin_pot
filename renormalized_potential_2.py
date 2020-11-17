# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 09:58:28 2019

@author: yadav
"""

#code to obtain normalized potential from f2(x)=V'(x) by numerical integration, see method 2 for non linear integral equation for f2(x) in notebook
# Note we add a constant of integration such that V_min=0 
import math
import cmath
import numpy
#import contour_integral
import integration_trapezoidal_nonuniformspacing
#import derivative_central_difference
#epsi = 1e-4
data1 = numpy.loadtxt("input_output_files/quadratic_V_rho=2.0/density/f_y_calculation_newton_rapson_method_2_CR_solvedproblem_gamma=0.5_rho=2.0_max_err_1e-10_18000points.txt",float)
#data2 = numpy.loadtxt("input_output_files//contour//nu2_contour_c1=1_c0=0.5_20000points_edit.txt",float)
#data3 = numpy.loadtxt("input_output_files//mapping//mapping_output_nu1_c1=1_c0=0.5_20000points.txt",float)
#contr = numpy.loadtxt("input_output_files//contour//nu_contour_c1=1_c0=0.5_20000points_edit.txt",float)
f_out=file("input_output_files/quadratic_V_rho=2.0/density/renormalized_potential_2_from_f2_y_method_2_CR_solvedproblem_gamma=0.5_max_err_1e-10_18000points.txt","w")
#x = data1[:, 0]    # Array Slicing: Accessing array's first column(i.e. index-zero column), https://jakevdp.github.io/PythonDataScienceHandbook/02.02-the-basics-of-numpy-arrays.html
#f2_x = data1[:, 1]    # Array Slicing: Accessing array's first column(i.e. index-zero column)
#deriv_V=f2_x/x    # since f2(x)=V'(x)x, see method 2 for non linear integral equation for f2(x) in notebook
V = numpy.zeros(len(data1),float) 

#for i in range(len(data1)-1):    # function len() on array gives no. of rows of array
for i in range(len(data1)):    # function len() on array gives no. of rows of array
#    V = integration_trapezoidal_nonuniformspacing.int_trap_non_uniform(x,deriv_V)
#    x1 = x[i]
#    f_out.write(str(x1)+" "+str(V)+'\n')
    j=len(data1)-1-i
    x = data1[j:, 0]    # Array Slicing: Accessing array's first column(i.e. index-zero column), https://jakevdp.github.io/PythonDataScienceHandbook/02.02-the-basics-of-numpy-arrays.html
    f2_x = data1[j:, 1]    # Array Slicing: Accessing array's first column(i.e. index-zero column)
    deriv_V = f2_x/x
    V[j] = (-1.0)*integration_trapezoidal_nonuniformspacing.int_trap_non_uniform(x,deriv_V)
#    V[j-1]=V[j]+(x[j-1]-x[j])*deriv_V[j]
#    x1 = x[j]
#    potential=V[j]
#    f_out.write(str(x1)+" "+str(potential)+'\n')
    print i
    
#print numpy.amin(V)  
V = V - numpy.amin(V)    # to add a constant of integration such that V_min=0 
    
for i in range(len(V)):
#    j=len(V)-1-i
    x1 = x[i]
    potential=V[i]
    f_out.write(str(x1)+" "+str(potential)+'\n')
        
f_out.close()    # () at the end is necessary to close the file     