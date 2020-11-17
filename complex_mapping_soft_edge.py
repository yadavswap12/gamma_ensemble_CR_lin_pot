# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 23:42:12 2019

@author: yadav
"""

# code to generate complex mapping J(s) for soft edge

#https://stackoverflow.com/questions/16669428/process-very-large-20gb-text-file-line-by-line
#modified from https://stackoverflow.com/questions/30216573/reading-specific-columns-from-a-text-file-in-python
#https://docs.python.org/3.3/tutorial/inputoutput.html
import math
import cmath
import numpy
theta = 2.0
rho = 2.0    # V(x)=rho*x 
#c = theta/rho    # eqn.(4.21) in C.R.

c1 = 0.747538234904    # renormalized c1 for gamma=0.8 computed from python file 'renormalized_joukowsky_parameters_c1_c0.py'
c0 = 0.853867428121    # renormalized c0 for gamma=0.8 computed from python file 'renormalized_joukowsky_parameters_c1_c0.py' 

#c0 = 1.0/2
f_in=file("input_output_files/V=rho_x/contour/nu1_contour_c1=0.5_c0=2.0_theta=2.0_soft_edge_18000points.txt","r")
f_out=file("input_output_files/V=rho_x/mapping/mapping_output_nu1_c1=0.5_c0=2.0_theta=2.0_soft_edge_18000points.txt","w")
lines=f_in.readlines()
#i=0
for data in lines:
    x= data.split(' ')[0]  
    y= data.split(' ')[1]  
    #z(i)=x(i)+y(i)j    gives syntax error
    z=complex(float(x),float(y))    
    x_nu1 = (c1*z+c0)*(((z+1.0)/z)**(1.0/theta))
    #f_out.write(str(x_nu1.real)+'\n')
    f_out.write(str(x_nu1.real)+" "+str(x_nu1.imag)+'\n')
    #i+=1
f_in.close()    # () at the end is necessary to close the file
f_out.close()    # () at the end is necessary to close the file

f_in=file("input_output_files/V=rho_x/contour/nu2_contour_c1=0.5_c0=2.0_theta=2.0_soft_edge_18000points.txt","r")
f_out=file("input_output_files/V=rho_x/mapping/mapping_output_nu2_c1=0.5_c0=2.0_theta=2.0_soft_edge_18000points.txt","w")
lines=f_in.readlines()
#i=0
for data in lines:
    x= data.split(' ')[0]  
    y= data.split(' ')[1]  
    #z(i)=x(i)+y(i)j    gives syntax error
    z=complex(float(x),float(y))
    x_nu2 = (c1*z+c0)*(((z+1.0)/z)**(1.0/theta))
    #f_out.write(str(x_nu2.real)+'\n')
    f_out.write(str(x_nu2.real)+" "+str(x_nu2.imag)+'\n')
    #i+=1
f_in.close()    # () at the end is necessary to close the file
f_out.close()    # () at the end is necessary to close the file