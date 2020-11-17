# -*- coding: utf-8 -*-
"""
Created on Tue Nov 27 10:38:00 2018

@author: yadav
"""

# code to generate nu1 mapping
import math
import cmath
import numpy
import matplotlib    #pylab is submodule in matplotlib
data1 = numpy.loadtxt("input_output_files/potential_linear_slope_2/mapping/mapping_output_nu1_gamma=0.8_theta=2.0_soft_edge_18000points_var_step_size.txt",float)
data2 = numpy.loadtxt("input_output_files/potential_linear_slope_2/mapping/mapping_output_nu2_gamma=0.8_theta=2.0_soft_edge_18000points_var_step_size.txt",float)
x1 = data1[:,0]
y1 = data1[:,1]
x2 = data2[:,0]
y2 = data2[:,1]
matplotlib.pylab.plot(x1,y1)
matplotlib.pylab.plot(x2,y2)
matplotlib.pylab.show()