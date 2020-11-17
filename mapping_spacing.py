# -*- coding: utf-8 -*-
"""
Created on Sun Jan 27 18:41:53 2019

@author: yadav
"""

#code for calculation of spacing between mapped points on real axis by joukowsky transformation.
#this code was needed to check if there was sharp difference(peak) in spacing of points as we have a sharp peak in density plot for soft edge joukowsky transformation from CR.   
 
import math
import cmath
import numpy
import contour_integral
import matplotlib    #pylab is submodule in matplotlib
#matplotlib.use('Agg')    # add this line to use code in linux

data1 = numpy.loadtxt("input_output_files/potential_linear_slope_2/mapping/mapping_output_nu1_gamma=0.8_theta=2.0_soft_edge_18000points_var_step_size.txt",float)
#data1 = numpy.loadtxt("input_output_files/potential_linear_slope_2/mapping/mapping_output_nu1_gamma=0.8_theta=2.0_18000points.txt",float)
f_out=file("input_output_files/potential_linear_slope_2/mapping/mapping_spacing_nu1_gamma=0.8_theta=2.0_soft_edge_18000points_var_step_size.txt","w")

for i in range(len(data1)-1):    # function len() on array gives no. of rows of array
    x = data1[i,0]
    spacing = data1[i+1,0] - data1[i,0]
    f_out.write(str(x)+" "+str(spacing)+'\n')
f_out.close()    # () at the end is necessary to close the file 

data2 = numpy.loadtxt("input_output_files/potential_linear_slope_2/mapping/mapping_spacing_nu1_gamma=0.8_theta=2.0_soft_edge_18000points_var_step_size.txt",float)
#data2 = numpy.loadtxt("input_output_files/potential_linear_slope_2/mapping/spacing_gamma=0.8_theta=2.0_18000points.txt",float)
#data3 = numpy.loadtxt("input_output_files/potential_linear_slope_2/density/renormalized_density_psi_method_2_epsi=1e-4_gamma=0.8_theta=2.0_rho=2.0_soft_edge_18000points.txt",float)
#data3 = numpy.loadtxt("input_output_files/potential_linear_slope_2/density/renormalized_density_psi_method_2_epsi=1e-4_gamma=0.8_theta=2.0_18000points_edited.txt",float)

x2 = data2[:,0]
y2 = data2[:,1]
#x3 = data3[:,0]
#y3 = data3[:,1]     

matplotlib.pylab.plot(x2,y2,"r")
#matplotlib.pylab.ylim(-0.0005,0)
matplotlib.pylab.show() 

#matplotlib.pylab.plot(x3,y3,"g")
#matplotlib.pylab.show() 
    
