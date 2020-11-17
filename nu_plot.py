# -*- coding: utf-8 -*-
"""
Created on Tue Nov 27 10:10:42 2018

@author: yadav
"""


# code to plot nu
import math
import cmath
import numpy
import matplotlib    #pylab is submodule in matplotlib
#import matplotlib.pylab 
#import pylab    # Use this step to run in anaconda
data = numpy.loadtxt("input_output_files/alpha_8.0_rho_2.0/gamma_0.99/density/phi_V_int_test_gamma=0.99_max_err_1e-10_theta=2.0_18000points_correction3.txt",float)
x = data[:,0]    # The plain [:] operator slices from beginning to end-1, Eg. a[:,0:n] gives all the rows from 0th column to n-1 column 
y = data[:,1]
#matplotlib.pylab.xlim(0.3,0.34)
matplotlib.pylab.plot(x,y)
matplotlib.pylab.ylim(-1,1)
#matplotlib.pyplot.plot(x,y)
#pylab.plot(x,y)    # Use this step to run in anaconda