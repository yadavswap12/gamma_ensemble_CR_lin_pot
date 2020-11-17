# -*- coding: utf-8 -*-
"""
Created on Tue Nov 27 10:08:32 2018

@author: yadav
"""

# code to generate contour nu, nu1 and nu2 for J(s)=c(s+1)((s+1)/s)**(1/theta)
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
#c = math.sqrt((1.0+theta)*(theta**2.0)/((2.0*rho)*(theta+2.0)))    # For V(x)=rho*x**2, see notes. 
 
#f_out=file("input_output_files/V=rho_x^2/contour/nu_contour_c=0.866_theta=2.0_rho=2.0_18000points.txt","w")
#f_out1=file("input_output_files/V=rho_x^2/contour/nu1_contour_c=0.866_theta=2.0_rho=2.0_18000points.txt","w")
#f_out2=file("input_output_files/V=rho_x^2/contour/nu2_contour_c=0.866_theta=2.0_rho=2.0_18000points.txt","w")

f_out=file("input_output_files/rho_2.0/gamma_"+str(gamma)+"/contour/nu_contour_c="+str(c_short)+"_theta=2.0_18000points_iter"+str(iteration)+".txt","w")
f_out1=file("input_output_files/rho_2.0/gamma_"+str(gamma)+"/contour/nu1_contour_c="+str(c_short)+"_theta=2.0_18000points_iter"+str(iteration)+".txt","w")
f_out2=file("input_output_files/rho_2.0/gamma_"+str(gamma)+"/contour/nu2_contour_c="+str(c_short)+"_theta=2.0_18000points_iter"+str(iteration)+".txt","w")


#for i in range(1,18000,1):    #only integer step argument is available for range() hence 0 to 18000. to get 180 degrees will later divide by 100
#    r = (math.tan((math.radians(i/100.0)/(1+theta))))/(math.sin(math.radians(i/100.0))-math.cos(math.radians(i/100.0))*math.tan((math.radians(i/100.0)/(1+theta))))    
#    x = r*math.cos(math.radians(i/100.0))
#    y = r*math.sin(math.radians(i/100.0))
#    #f_out.write(str(r)+" "+str(x)+" "+str(y)+" "+str(i)+'\n')
#    f_out.write(str(x)+" "+str(y)+'\n')
#    f_out1.write(str(x)+" "+str(y)+'\n')
##f_out.close()    # () at the end is necessary to close the file
#f_out1.close()    # () at the end is necessary to close the file
#for i in range(18000,1,-1):    #only integer step argument is available for range() hence 0 to 18000. to get 180 degrees will later divide by 100
#    r = (math.tan((math.radians(i/100.0)/(1+theta))))/(math.sin(math.radians(i/100.0))-math.cos(math.radians(i/100.0))*math.tan((math.radians(i/100.0)/(1+theta))))    
#    x = r*math.cos(math.radians(i/100.0))
#    y = -r*math.sin(math.radians(i/100.0))
#    #f_out.write(str(r)+" "+str(x)+" "+str(y)+" "+str(i)+'\n')
#    f_out.write(str(x)+" "+str(y)+'\n')
#    f_out2.write(str(x)+" "+str(y)+'\n')
#f_out.close()    # () at the end is necessary to close the file    
#f_out2.close()    # () at the end is necessary to close the file

for i in range(1,18000,1):    #only integer step argument is available for range() hence 0 to 36000. to get 180 degrees will later divide by 200
    r = (math.tan((math.radians(i/100.0)/(1.0+theta))))/(math.sin(math.radians(i/100.0))-math.cos(math.radians(i/100.0))*math.tan((math.radians(i/100.0)/(1.0+theta))))    
    x = r*math.cos(math.radians(i/100.0))
    y = r*math.sin(math.radians(i/100.0))
    #f_out.write(str(r)+" "+str(x)+" "+str(y)+" "+str(i)+'\n')
    f_out.write(str(x)+" "+str(y)+'\n')
    f_out1.write(str(x)+" "+str(y)+'\n')
#f_out.close()    # () at the end is necessary to close the file
f_out1.close()    # () at the end is necessary to close the file
#for i in range(18000,1,-1):    #only integer step argument is available for range() hence 0 to 36000. to get 180 degrees will later divide by 200
for i in range(17999,0,-1):    #only integer step argument is available for range() hence 0 to 36000. to get 180 degrees will later divide by 200
    r = (math.tan((math.radians(i/100.0)/(1.0+theta))))/(math.sin(math.radians(i/100.0))-math.cos(math.radians(i/100.0))*math.tan((math.radians(i/100.0)/(1.0+theta))))    
    x = r*math.cos(math.radians(i/100.0))
    y = -r*math.sin(math.radians(i/100.0))
    #f_out.write(str(r)+" "+str(x)+" "+str(y)+" "+str(i)+'\n')
    f_out.write(str(x)+" "+str(y)+'\n')
    f_out2.write(str(x)+" "+str(y)+'\n')
f_out.close()    # () at the end is necessary to close the file    
f_out2.close()    # () at the end is necessary to close the file          
   