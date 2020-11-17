# -*- coding: utf-8 -*-
"""
Created on Tue Nov 27 10:12:39 2018

@author: yadav
"""

# code to generate complex mapping J(s)

#https://stackoverflow.com/questions/16669428/process-very-large-20gb-text-file-line-by-line
#modified from https://stackoverflow.com/questions/30216573/reading-specific-columns-from-a-text-file-in-python
#https://docs.python.org/3.3/tutorial/inputoutput.html
import math
import cmath
import numpy

theta = 2.0
rho = 2.0    # V(x)=rho*x+alpha*x**(1.0/alpha), alpha>1
#alpha = 8.0    # V(x)=rho*x+alpha*x**(1.0/alpha), alpha>1
c = 0.4047464827825728    # from self-consistent calculation using jouwkowsky_parameters_c_selfconsistent_calculation_CR_problem.py for alpha-log(x) potential.
c_short = 0.404
gamma = 0.1
iteration = 3.0

#theta = 2.0    # parameters from solved C.R. problem
#rho = 2.0    # V(x)=rho*x**2
#c = math.sqrt((1.0+theta)*(theta**2)/((2*rho)*(theta+2)))    # For V(x)=rho*x**2, see notes.

#c0 = 1.0/2

#f_in=file("input_output_files/V=rho_x^2/contour/nu1_contour_c=0.866_theta=2.0_rho=2.0_18000points.txt","r")
#f_out=file("input_output_files/V=rho_x^2/mapping/mapping_output_nu1_c=0.866_theta=2.0_rho=2.0_18000points.txt","w")

f_in=file("input_output_files/rho_2.0/gamma_"+str(gamma)+"/contour/nu1_contour_c="+str(c_short)+"_theta=2.0_18000points_iter"+str(iteration)+".txt","r")
f_out=file("input_output_files/rho_2.0/gamma_"+str(gamma)+"/mapping/mapping_output_nu1_c="+str(c_short)+"_theta=2.0_18000points_iter"+str(iteration)+".txt","w")


lines=f_in.readlines()
#i=0
for data in lines:
    x= data.split(' ')[0]  
    y= data.split(' ')[1]  
    #z(i)=x(i)+y(i)j    #gives syntax error, use z(i)=x(i)+y(i)*1j 
    z=complex(float(x),float(y))    
    x_nu1 = c*(z+1)*(((z+1)/z)**(1/theta))
    #f_out.write(str(x_nu1.real)+'\n')
    f_out.write(str(x_nu1.real)+" "+str(x_nu1.imag)+'\n')
    #i+=1
f_in.close()    # () at the end is necessary to close the file
f_out.close()    # () at the end is necessary to close the file

#f_in=file("input_output_files/V=rho_x^2/contour/nu2_contour_c=0.866_theta=2.0_rho=2.0_18000points.txt","r")
#f_out=file("input_output_files/V=rho_x^2/mapping/mapping_output_nu2_c=0.866_theta=2.0_rho=2.0_18000points.txt","w")

f_in=file("input_output_files/rho_2.0/gamma_"+str(gamma)+"/contour/nu2_contour_c="+str(c_short)+"_theta=2.0_18000points_iter"+str(iteration)+".txt","r")
f_out=file("input_output_files/rho_2.0/gamma_"+str(gamma)+"/mapping/mapping_output_nu2_c="+str(c_short)+"_theta=2.0_18000points_iter"+str(iteration)+".txt","w")
lines=f_in.readlines()
#i=0
for data in lines:
    x= data.split(' ')[0]  
    y= data.split(' ')[1]  
    #z(i)=x(i)+y(i)j    gives syntax error
    z=complex(float(x),float(y))
    x_nu2 = c*(z+1)*(((z+1)/z)**(1/theta))
    #f_out.write(str(x_nu2.real)+'\n')
    f_out.write(str(x_nu2.real)+" "+str(x_nu2.imag)+'\n')
    #i+=1
f_in.close()    # () at the end is necessary to close the file
f_out.close()    # () at the end is necessary to close the file