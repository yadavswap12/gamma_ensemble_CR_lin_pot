# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 10:16:58 2019

@author: yadav
"""

#code for calculation of spacing between contour points for joukowsky transformation.
#this code was needed to check if there was sharp difference(peak) in spacing of points as we have a sharp peak in density plot for soft edge joukowsky transformation from CR.   
 
import math
import cmath
import numpy
import contour_integral
import matplotlib    #pylab is submodule in matplotlib
#matplotlib.use('Agg')    # add this line to use code in linux

data1 = numpy.loadtxt("input_output_files/rho_2.0/contour/nu_contour_c=1_theta=2.0_18000points.txt",float)
#data1 = numpy.loadtxt("input_output_files/potential_linear_slope_2/mapping/mapping_output_nu1_gamma=0.8_theta=2.0_18000points.txt",float)
f_out=file("input_output_files/rho_2.0/contour/contour_spacing_nu_contour_c=1_theta=2.0_18000points.txt","w")

for i in range(len(data1)-1):    # function len() on array gives no. of rows of array
    x = data1[i,0]
    y = data1[i,1]
    phi = i
    x_spacing = data1[i+1,0] - data1[i,0]
    y_spacing = data1[i+1,1] - data1[i,1]    
    f_out.write(str(phi)+" "+str(x_spacing)+" "+str(y_spacing)+'\n')
f_out.close()    # () at the end is necessary to close the file 

data2 = numpy.loadtxt("input_output_files/rho_2.0/contour/contour_spacing_nu_contour_c=1_theta=2.0_18000points.txt",float)
#data2 = numpy.loadtxt("input_output_files/potential_linear_slope_2/mapping/spacing_gamma=0.8_theta=2.0_18000points.txt",float)
#data3 = numpy.loadtxt("input_output_files/potential_linear_slope_2/density/renormalized_density_psi_method_2_epsi=1e-4_gamma=0.8_theta=2.0_rho=2.0_soft_edge_18000points.txt",float)
#data3 = numpy.loadtxt("input_output_files/potential_linear_slope_2/density/renormalized_density_psi_method_2_epsi=1e-4_gamma=0.8_theta=2.0_18000points_edited.txt",float)

x2 = data2[:,0]
y2 = data2[:,1]
y3 = data2[:,2]
#x3 = data3[:,0]
#y3 = data3[:,1] 

    

matplotlib.pylab.plot(x2,y2,"r")
matplotlib.pylab.ylabel('x-spacing')
matplotlib.pylab.xlabel('$\phi$')
#matplotlib.pylab.ylim(-0.005,0.005)
#matplotlib.pylab.xlim(94.6,95)
#matplotlib.pylab.xlim(90,95)
matplotlib.pylab.show()

matplotlib.pylab.plot(x2,y3,"g")
matplotlib.pylab.ylabel('y-spacing')
matplotlib.pylab.xlabel('$\phi$')
#matplotlib.pylab.ylim(-0.005,0)
#matplotlib.pylab.xlim(94.6,95)
#matplotlib.pylab.xlim(90,95)
matplotlib.pylab.show()  

#matplotlib.pylab.plot(x3,y3,"g")
#matplotlib.pylab.show() 