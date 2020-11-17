# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 17:56:26 2019

@author: yadav
"""
# code to generate contour nu, nu1 and nu2 for J(s)=(c1s+c0)((s+1)/s)**(1/theta) with variable step-size using limit of normal distribution as dirac-delta function. See https://en.wikipedia.org/wiki/Dirac_delta_function
# Note that we choose to increase or decrease step-size depending on output of contour_spacing.py 
import math
import cmath
import numpy
theta = 2.0
rho = 2.0    # V(x)=rho*x
#c = theta/rho    # eqn.(4.21) in C.R.
c1 = 0.873607034265    # renormalized c1 for gamma=0.9 computed from python file 'renormalized_joukowsky_parameters_c1_c0.py'
c0 = 0.91688096976    # renormalized c0 for gamma=0.9 computed from python file 'renormalized_joukowsky_parameters_c1_c0.py' 

#aa0=800000.0
#aa=50000.0
x0=10502
#aa0=80000.0
#aa=5000.0

#aa0=4000000.0
#aa=100000.0

aa0=40000000.0
aa=1000000.0

f_out=file("input_output_files/potential_linear_slope_2/contour/nu_contour_gamma=0.9_theta=2.0_soft_edge_18000points_var_step_size2.txt","w")
f_out1=file("input_output_files/potential_linear_slope_2/contour/nu1_contour_gamma=0.9_theta=2.0_soft_edge_18000points_var_step_size2.txt","w")
f_out2=file("input_output_files/potential_linear_slope_2/contour/nu2_contour_gamma=0.9_theta=2.0_soft_edge_18000points_var_step_size2.txt","w")

for i in range(1,18000,1):    #only integer step argument is available for range() hence 0 to 18000. to get 180 degrees will later divide by 100

#    ii=i*(1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((i-x0)**2.0)/aa))

    for j in range(i*int((1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((i-x0)**2.0)/aa))),(i+1)*int((1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((i-x0)**2.0)/aa))),1):    # to create smaller step-size.    #only integer step argument is available for range() hence 0 to 18000. to get 180 degrees will later divide by 100
        cos_phi =  math.cos(math.radians(j/((100.0)*int((1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((i-x0)**2.0)/aa))))))
        sin_phi =  math.sin(math.radians(j/((100.0)*int((1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((i-x0)**2.0)/aa))))))
        tan_phi =  math.tan(math.radians(j/((100.0)*int((1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((i-x0)**2.0)/aa))))))
        tan_phibytheta =  math.tan(math.radians(j/((100.0)*int((1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((i-x0)**2.0)/aa))))/theta))
        tan_minusphibytheta =  math.tan(math.radians(-j/((100.0)*int((1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((i-x0)**2.0)/aa)))))/theta)
#        phi = math.atan((-sin_phi_prime)/())
        a = 1.0
#        b = ((c1/c0)*(cos_phi*tan_phibytheta-sin_phi)*tan_phibytheta-2*cos_phi)
#        b = (c1/c0)*(cos_phi-(sin_phi)/tan_phibytheta)-2.0*cos_phi
        b = ((c1*sin_phi)/(c0*tan_minusphibytheta) + (c1*cos_phi)/(c0) - 2.0*cos_phi)
#        c = 1-(c1/c0)*(tan_phibytheta**2)
        c = 1.0-(c1/c0)
#        print i
        if(b**2.0-4.0*a*c>=0):
            r = (-b+(b**2.0-4.0*a*c)**(1.0/2))/(2.0*a)
            x = (r*cos_phi-1)/(r**2.0-2.0*r*cos_phi+1)
            y = (r*sin_phi)/(r**2.0-2.0*r*cos_phi+1)
            f_out.write(str(x)+" "+str(y)+" "+str(j/((100.0)*int((1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((i-x0)**2.0)/aa)))))+'\n')
            f_out1.write(str(x)+" "+str(y)+" "+str(j/((100.0)*int((1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((i-x0)**2.0)/aa)))))+'\n')
        
        else:
            break

#    print('i', i)
        
    if(b**2.0-4.0*a*c>=0):
        continue
    else:
        break

#print i,j            

for k in range(i-1,0,-1):    #only integer step argument is available for range() hence 0 to 18000. to get 180 degrees will later divide by 100

    for j in range(k*int((1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((k-x0)**2.0)/aa))),(k-1)*int((1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((k-x0)**2.0)/aa))),-1):    # to create smaller step-size.    #only integer step argument is available for range() hence 0 to 18000. to get 180 degrees will later divide by 100
        cos_phi =  math.cos(math.radians(j/((100.0)*int((1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((k-x0)**2.0)/aa))))))
        sin_phi =  math.sin(math.radians(j/((100.0)*int((1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((k-x0)**2.0)/aa))))))
        tan_phi =  math.tan(math.radians(j/((100.0)*int((1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((k-x0)**2.0)/aa))))))
        tan_phibytheta =  math.tan(math.radians(j/((100.0)*int((1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((k-x0)**2.0)/aa)))))/theta)
        tan_minusphibytheta =  math.tan(math.radians(-j/((100.0)*int((1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((k-x0)**2.0)/aa)))))/theta)
#        phi = math.atan((-sin_phi_prime)/())
        a = 1.0
#        b = ((c1/c0)*(cos_phi*tan_phibytheta-sin_phi)*tan_phibytheta-2*cos_phi)
#        b = (c1/c0)*(cos_phi-(sin_phi)/tan_phibytheta)-2.0*cos_phi
        b = ((c1*sin_phi)/(c0*tan_minusphibytheta) + (c1*cos_phi)/(c0) - 2.0*cos_phi)
#        c = 1-(c1/c0)*(tan_phibytheta**2)
        c = 1.0-(c1/c0)
#        print i
        if(b**2.0-4.0*a*c>=0):
            r = (-b-(b**2.0-4.0*a*c)**(1.0/2))/(2.0*a)
            x = (r*cos_phi-1)/(r**2.0-2.0*r*cos_phi+1)
            y = (r*sin_phi)/(r**2.0-2.0*r*cos_phi+1)
            f_out.write(str(x)+" "+str(y)+" "+str(j/((100.0)*int((1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((k-x0)**2.0)/aa)))))+'\n')
            f_out1.write(str(x)+" "+str(y)+" "+str(j/((100.0)*int((1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((k-x0)**2.0)/aa)))))+'\n')
        
        else:
            break

#    print('k', k)    
    
    if(b**2.0-4.0*a*c>=0):
        continue
    else:
        break
      
#f_out.close()    # () at the end is necessary to close the file
f_out1.close()    # () at the end is necessary to close the file


for i in range(1,18000,1):    #only integer step argument is available for range() hence 0 to 18000. to get 180 degrees will later divide by 100

    for j in range(i*int((1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((i-x0)**2.0)/aa))),(i+1)*int((1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((i-x0)**2.0)/aa))),1):    # to create smaller step-size.    #only integer step argument is available for range() hence 0 to 18000. to get 180 degrees will later divide by 100
        cos_phi =  math.cos(math.radians(j/((100.0)*int((1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((i-x0)**2.0)/aa))))))
        sin_phi =  math.sin(math.radians(j/((100.0)*int((1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((i-x0)**2.0)/aa))))))
        tan_phi =  math.tan(math.radians(j/((100.0)*int((1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((i-x0)**2.0)/aa))))))
        tan_phibytheta =  math.tan(math.radians(j/((100.0)*int((1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((i-x0)**2.0)/aa))))/theta))
        tan_minusphibytheta =  math.tan(math.radians(-j/((100.0)*int((1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((i-x0)**2.0)/aa)))))/theta)
#        phi = math.atan((-sin_phi_prime)/())
        a = 1.0
#        b = ((c1/c0)*(cos_phi*tan_phibytheta-sin_phi)*tan_phibytheta-2*cos_phi)
#        b = (c1/c0)*(cos_phi-(sin_phi)/tan_phibytheta)-2.0*cos_phi
        b = ((c1*sin_phi)/(c0*tan_minusphibytheta) + (c1*cos_phi)/(c0) - 2.0*cos_phi)
#        c = 1-(c1/c0)*(tan_phibytheta**2)
        c = 1.0-(c1/c0)
#        print i
        if(b**2.0-4.0*a*c>=0):
            r = (-b-(b**2.0-4.0*a*c)**(1.0/2))/(2.0*a)
            x = (r*cos_phi-1)/(r**2.0-2.0*r*cos_phi+1)
            y = (-r*sin_phi)/(r**2.0-2.0*r*cos_phi+1)
            f_out.write(str(x)+" "+str(y)+" "+str(-j/((100.0)*int((1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((i-x0)**2.0)/aa)))))+'\n')
            f_out2.write(str(x)+" "+str(y)+" "+str(-j/((100.0)*int((1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((i-x0)**2.0)/aa)))))+'\n')
        
        else:
            break
        
#    print('i', i)        
        
    if(b**2.0-4.0*a*c>=0):
        continue
    else:
        break

for k in range(i-1,0,-1):    #only integer step argument is available for range() hence 0 to 18000. to get 180 degrees will later divide by 100

    for j in range(k*int((1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((k-x0)**2.0)/aa))),(k-1)*int((1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((k-x0)**2.0)/aa))),-1):    # to create smaller step-size.    #only integer step argument is available for range() hence 0 to 18000. to get 180 degrees will later divide by 100
        cos_phi =  math.cos(math.radians(j/((100.0)*int((1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((k-x0)**2.0)/aa))))))
        sin_phi =  math.sin(math.radians(j/((100.0)*int((1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((k-x0)**2.0)/aa))))))
        tan_phi =  math.tan(math.radians(j/((100.0)*int((1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((k-x0)**2.0)/aa))))))
        tan_phibytheta =  math.tan(math.radians(j/((100.0)*int((1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((k-x0)**2.0)/aa)))))/theta)
        tan_minusphibytheta =  math.tan(math.radians(-j/((100.0)*int((1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((k-x0)**2.0)/aa)))))/theta)
#        phi = math.atan((-sin_phi_prime)/())
        a = 1.0
#        b = ((c1/c0)*(cos_phi*tan_phibytheta-sin_phi)*tan_phibytheta-2*cos_phi)
#        b = (c1/c0)*(cos_phi-(sin_phi)/tan_phibytheta)-2.0*cos_phi
        b = ((c1*sin_phi)/(c0*tan_minusphibytheta) + (c1*cos_phi)/(c0) - 2.0*cos_phi)
#        c = 1-(c1/c0)*(tan_phibytheta**2)
        c = 1.0-(c1/c0)
#        print i
        if(b**2.0-4.0*a*c>=0):
            r = (-b+(b**2.0-4.0*a*c)**(1.0/2))/(2.0*a)
            x = (r*cos_phi-1)/(r**2.0-2.0*r*cos_phi+1)
            y = (-r*sin_phi)/(r**2.0-2.0*r*cos_phi+1)
            f_out.write(str(x)+" "+str(y)+" "+str(-j/((100.0)*int((1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((k-x0)**2.0)/aa)))))+'\n')
            f_out2.write(str(x)+" "+str(y)+" "+str(-j/((100.0)*int((1.0+(1.0/(aa*math.sqrt(math.pi)))*aa0*math.exp(-((k-x0)**2.0)/aa)))))+'\n')
        
        else:
            break
        
#    print('k', k)    
        
    if(b**2.0-4.0*a*c>=0):
        continue
    else:
        break
        
f_out.close()    # () at the end is necessary to close the file
f_out2.close()    # () at the end is necessary to close the file