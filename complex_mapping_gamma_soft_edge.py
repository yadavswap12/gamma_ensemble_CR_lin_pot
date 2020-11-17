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
alpha=8.0
rho = 2.0    # V(x)=gamma*(rho*x+alpha*x**(1.0/alpha)), alpha>1
gamma = 0.76
iteration = 11.0

#c = theta/rho    # eqn.(4.21) in C.R.
c1 = 0.75651787191897635    # renormalized c1 for gamma=0.5 computed from python file 'renormalized_joukowsky_parameters_c1_c0.py'
c0 = 0.76780260554911939    # renormalized c0 for gamma=0.5 computed from python file 'renormalized_joukowsky_parameters_c1_c0.py' 

#c0 = 1.0/2
f_in=file("input_output_files/alpha_"+str(alpha)+"_rho_"+str(rho)+"/gamma_"+str(gamma)+"/contour/nu1_contour_gamma="+str(gamma)+"_theta="+str(theta)+"_soft_edge_18000points_corrected3_iter"+str(iteration)+".txt","r")
f_out=file("input_output_files/alpha_"+str(alpha)+"_rho_"+str(rho)+"/gamma_"+str(gamma)+"/mapping/mapping_output_nu1_gamma="+str(gamma)+"_theta="+str(theta)+"_soft_edge_18000points_corrected3_iter"+str(iteration)+".txt","w")
lines=f_in.readlines()
#i=0
for data in lines:
    x= data.split(' ')[0]  
    y= data.split(' ')[1]  
    #z(i)=x(i)+y(i)j    gives syntax error
    z=complex(float(x),float(y))    
    x_nu1 = (c1*z+c0)*(((z+1)/z)**(1/theta))
    #f_out.write(str(x_nu1.real)+'\n')
    f_out.write(str(x_nu1.real)+" "+str(x_nu1.imag)+'\n')
    #i+=1
f_in.close()    # () at the end is necessary to close the file
f_out.close()    # () at the end is necessary to close the file

f_in=file("input_output_files/alpha_"+str(alpha)+"_rho_"+str(rho)+"/gamma_"+str(gamma)+"/contour/nu2_contour_gamma="+str(gamma)+"_theta="+str(theta)+"_soft_edge_18000points_corrected3_iter"+str(iteration)+".txt","r")
f_out=file("input_output_files/alpha_"+str(alpha)+"_rho_"+str(rho)+"/gamma_"+str(gamma)+"/mapping/mapping_output_nu2_gamma="+str(gamma)+"_theta="+str(theta)+"_soft_edge_18000points_corrected3_iter"+str(iteration)+".txt","w")
lines=f_in.readlines()
#i=0
for data in lines:
    x= data.split(' ')[0]  
    y= data.split(' ')[1]  
    #z(i)=x(i)+y(i)j    gives syntax error
    z=complex(float(x),float(y))
    x_nu2 = (c1*z+c0)*(((z+1)/z)**(1/theta))
    #f_out.write(str(x_nu2.real)+'\n')
    f_out.write(str(x_nu2.real)+" "+str(x_nu2.imag)+'\n')
    #i+=1
f_in.close()    # () at the end is necessary to close the file
f_out.close()    # () at the end is necessary to close the file