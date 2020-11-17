# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 11:21:21 2019

@author: yadav
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 09:14:32 2018

@author: yadav
"""

# code to generate contour nu, nu1 and nu2 for J(s)=(c1s+c0)((s+1)/s)**(1/theta) with variable step-size.
# Note that we choose to increase or decrease step-size depending on output of contour_spacing.py 
import math
import cmath
import numpy
theta = 2.0
rho = 2.0    # V(x)=rho*x
#c = theta/rho    # eqn.(4.21) in C.R.
c1 = 0.747538234904    # renormalized c1 for gamma=0.8 computed from python file 'renormalized_joukowsky_parameters_c1_c0.py'
c0 = 0.853867428121    # renormalized c0 for gamma=0.8 computed from python file 'renormalized_joukowsky_parameters_c1_c0.py' 

f_out=file("input_output_files/potential_linear_slope_2/contour/nu_contour_gamma=0.8_theta=2.0_soft_edge_18000points_var_step_size.txt","w")
f_out1=file("input_output_files/potential_linear_slope_2/contour/nu1_contour_gamma=0.8_theta=2.0_soft_edge_18000points_var_step_size.txt","w")
f_out2=file("input_output_files/potential_linear_slope_2/contour/nu2_contour_gamma=0.8_theta=2.0_soft_edge_18000points_var_step_size.txt","w")

for i in range(1,18000,1):    #only integer step argument is available for range() hence 0 to 18000. to get 180 degrees will later divide by 100

    if(i<9350):
        cos_phi =  math.cos(math.radians(i/100.0))
        sin_phi =  math.sin(math.radians(i/100.0))
        tan_phi =  math.tan(math.radians(i/100.0))
        tan_phibytheta =  math.tan(math.radians(i/100.0)/theta)
        tan_minusphibytheta =  math.tan(math.radians(-i/100.0)/theta)
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
            f_out.write(str(x)+" "+str(y)+" "+str(i/100.0)+'\n')
            f_out1.write(str(x)+" "+str(y)+" "+str(i/100.0)+'\n')
        
        else:
            break
    else:
        for j in range(i*10,(i+1)*10,1):    # to create smaller step-size.    #only integer step argument is available for range() hence 0 to 18000. to get 180 degrees will later divide by 100
            cos_phi =  math.cos(math.radians(j/1000.0))
            sin_phi =  math.sin(math.radians(j/1000.0))
            tan_phi =  math.tan(math.radians(j/1000.0))
            tan_phibytheta =  math.tan(math.radians(j/1000.0)/theta)
            tan_minusphibytheta =  math.tan(math.radians(-j/1000.0)/theta)
#            phi = math.atan((-sin_phi_prime)/())
            a = 1.0
#            b = ((c1/c0)*(cos_phi*tan_phibytheta-sin_phi)*tan_phibytheta-2*cos_phi)
#            b = (c1/c0)*(cos_phi-(sin_phi)/tan_phibytheta)-2.0*cos_phi
            b = ((c1*sin_phi)/(c0*tan_minusphibytheta) + (c1*cos_phi)/(c0) - 2.0*cos_phi)
#            c = 1-(c1/c0)*(tan_phibytheta**2)
            c = 1.0-(c1/c0)
#            print i
            if(b**2.0-4.0*a*c>=0):
                r = (-b+(b**2.0-4.0*a*c)**(1.0/2))/(2.0*a)
                x = (r*cos_phi-1)/(r**2.0-2.0*r*cos_phi+1)
                y = (r*sin_phi)/(r**2.0-2.0*r*cos_phi+1)
                f_out.write(str(x)+" "+str(y)+" "+str(j/1000.0)+'\n')
                f_out1.write(str(x)+" "+str(y)+" "+str(j/1000.0)+'\n')
        
            else:
                break

        if(b**2.0-4.0*a*c>=0):
            continue
        else:
            break

#        else:
#            break       

print i,j        
    
#for k in range(i-1,0,-1):    #only integer step argument is available for range() hence 0 to 18000. to get 180 degrees will later divide by 100
for k in range(i,0,-1):    #only integer step argument is available for range() hence 0 to 18000. to get 180 degrees will later divide by 100

    if(k<9350):
        cos_phi =  math.cos(math.radians(k/100.0))
        sin_phi =  math.sin(math.radians(k/100.0))
        tan_phi =  math.tan(math.radians(k/100.0))
        tan_phibytheta =  math.tan(math.radians(k/100.0)/theta)
        tan_minusphibytheta =  math.tan(math.radians(-k/100.0)/theta)
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
            f_out.write(str(x)+" "+str(y)+" "+str(k/100.0)+'\n')
            f_out1.write(str(x)+" "+str(y)+" "+str(k/100.0)+'\n')
        
        else:
            break
    else:
        for j in range(k*10,(k-1)*10,-1):    # to create smaller step-size.    #only integer step argument is available for range() hence 0 to 18000. to get 180 degrees will later divide by 100
            cos_phi =  math.cos(math.radians(j/1000.0))
            sin_phi =  math.sin(math.radians(j/1000.0))
            tan_phi =  math.tan(math.radians(j/1000.0))
            tan_phibytheta =  math.tan(math.radians(j/1000.0)/theta)
            tan_minusphibytheta =  math.tan(math.radians(-j/1000.0)/theta)
#            phi = math.atan((-sin_phi_prime)/())
            a = 1.0
#            b = ((c1/c0)*(cos_phi*tan_phibytheta-sin_phi)*tan_phibytheta-2*cos_phi)
#            b = (c1/c0)*(cos_phi-(sin_phi)/tan_phibytheta)-2.0*cos_phi
            b = ((c1*sin_phi)/(c0*tan_minusphibytheta) + (c1*cos_phi)/(c0) - 2.0*cos_phi)
#            c = 1-(c1/c0)*(tan_phibytheta**2)
            c = 1.0-(c1/c0)
#            print i
            if(b**2.0-4.0*a*c>=0):
                r = (-b-(b**2.0-4.0*a*c)**(1.0/2))/(2.0*a)
                x = (r*cos_phi-1)/(r**2.0-2.0*r*cos_phi+1)
                y = (r*sin_phi)/(r**2.0-2.0*r*cos_phi+1)
                f_out.write(str(x)+" "+str(y)+" "+str(j/1000.0)+'\n')
                f_out1.write(str(x)+" "+str(y)+" "+str(j/1000.0)+'\n')
        
            else:
                break

        if(b**2.0-4.0*a*c>=0):
            continue
        else:
            break
            
#        else:
#            break
    
        
#f_out.close()    # () at the end is necessary to close the file
f_out1.close()    # () at the end is necessary to close the file


for i in range(1,18000,1):    #only integer step argument is available for range() hence 0 to 18000. to get 180 degrees will later divide by 100

    if(i<9350):
        cos_phi =  math.cos(math.radians(i/100.0))
        sin_phi =  math.sin(math.radians(i/100.0))
        tan_phi =  math.tan(math.radians(i/100.0))
        tan_phibytheta =  math.tan(math.radians(i/100.0)/theta)
        tan_minusphibytheta =  math.tan(math.radians(-i/100.0)/theta)
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
            f_out.write(str(x)+" "+str(y)+" "+str(-i/100.0)+'\n')
            f_out2.write(str(x)+" "+str(y)+" "+str(-i/100.0)+'\n')
        
        else:
            break

    else:
        for j in range(i*10,(i+1)*10,1):    # to create smaller step-size.    #only integer step argument is available for range() hence 0 to 18000. to get 180 degrees will later divide by 100
            cos_phi =  math.cos(math.radians(j/1000.0))
            sin_phi =  math.sin(math.radians(j/1000.0))
            tan_phi =  math.tan(math.radians(j/1000.0))
            tan_phibytheta =  math.tan(math.radians(j/1000.0)/theta)
            tan_minusphibytheta =  math.tan(math.radians(-j/1000.0)/theta)
#            phi = math.atan((-sin_phi_prime)/())
            a = 1.0
#            b = ((c1/c0)*(cos_phi*tan_phibytheta-sin_phi)*tan_phibytheta-2*cos_phi)
#            b = (c1/c0)*(cos_phi-(sin_phi)/tan_phibytheta)-2.0*cos_phi
            b = ((c1*sin_phi)/(c0*tan_minusphibytheta) + (c1*cos_phi)/(c0) - 2.0*cos_phi)
#            c = 1-(c1/c0)*(tan_phibytheta**2)
            c = 1.0-(c1/c0)
#            print i
            if(b**2.0-4.0*a*c>=0):
                r = (-b-(b**2.0-4.0*a*c)**(1.0/2))/(2.0*a)
                x = (r*cos_phi-1)/(r**2.0-2.0*r*cos_phi+1)
                y = (-r*sin_phi)/(r**2.0-2.0*r*cos_phi+1)
                f_out.write(str(x)+" "+str(y)+" "+str(j/1000.0)+'\n')
                f_out2.write(str(x)+" "+str(y)+" "+str(j/1000.0)+'\n')
        
            else:
                break

        if(b**2.0-4.0*a*c>=0):
            continue
        else:
            break        
    
#for k in range(i-1,0,-1):    #only integer step argument is available for range() hence 0 to 18000. to get 180 degrees will later divide by 100
for k in range(i,0,-1):    #only integer step argument is available for range() hence 0 to 18000. to get 180 degrees will later divide by 100

    if(k<9350):
        cos_phi =  math.cos(math.radians(k/100.0))
        sin_phi =  math.sin(math.radians(k/100.0))
        tan_phi =  math.tan(math.radians(k/100.0))
        tan_phibytheta =  math.tan(math.radians(k/100.0)/theta)
        tan_minusphibytheta =  math.tan(math.radians(-k/100.0)/theta)
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
            f_out.write(str(x)+" "+str(y)+" "+str(-k/100.0)+'\n')
            f_out2.write(str(x)+" "+str(y)+" "+str(-k/100.0)+'\n')
        
        else:
            break
        
    else:
        for j in range(k*10,(k-1)*10,-1):    # to create smaller step-size.    #only integer step argument is available for range() hence 0 to 18000. to get 180 degrees will later divide by 100
            cos_phi =  math.cos(math.radians(j/1000.0))
            sin_phi =  math.sin(math.radians(j/1000.0))
            tan_phi =  math.tan(math.radians(j/1000.0))
            tan_phibytheta =  math.tan(math.radians(j/1000.0)/theta)
            tan_minusphibytheta =  math.tan(math.radians(-j/1000.0)/theta)
#            phi = math.atan((-sin_phi_prime)/())
            a = 1.0
#            b = ((c1/c0)*(cos_phi*tan_phibytheta-sin_phi)*tan_phibytheta-2*cos_phi)
#            b = (c1/c0)*(cos_phi-(sin_phi)/tan_phibytheta)-2.0*cos_phi
            b = ((c1*sin_phi)/(c0*tan_minusphibytheta) + (c1*cos_phi)/(c0) - 2.0*cos_phi)
#            c = 1-(c1/c0)*(tan_phibytheta**2)
            c = 1.0-(c1/c0)
#            print i
            if(b**2.0-4.0*a*c>=0):
                r = (-b+(b**2.0-4.0*a*c)**(1.0/2))/(2.0*a)
                x = (r*cos_phi-1)/(r**2.0-2.0*r*cos_phi+1)
                y = (-r*sin_phi)/(r**2.0-2.0*r*cos_phi+1)
                f_out.write(str(x)+" "+str(y)+" "+str(j/1000.0)+'\n')
                f_out2.write(str(x)+" "+str(y)+" "+str(j/1000.0)+'\n')
        
            else:
                break

        if(b**2.0-4.0*a*c>=0):
            continue
        else:
            break
            
f_out.close()    # () at the end is necessary to close the file
f_out2.close()    # () at the end is necessary to close the file