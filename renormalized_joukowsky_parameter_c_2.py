# -*- coding: utf-8 -*-
"""
Created on Mon Jan 07 16:50:40 2019

@author: yadav
"""

#code for computation of parameter 'c' of Joukowsky transform self consistently for application of NUS paper method to solved example in CR paper. See notebook for details
import math
import cmath
import numpy
import contour_integral

b0=0.31415    # parameters of polynomial fit for f1(x), f1(x)=b0 + b1x + b2x^2 + b3x^3 + b4x^4 
b1=8.89297    # gamma = 0.75
b2=-35.76914
b3=149.59459
b4=-369.06716
b5=541.69376
b6=-464.46166
b7=214.55806
b8=-41.21924

c0=0.58341
c1=3.95828

#gamma=0.77 c=0.75964909760813049    #corrected3 iteration11
#gamma=0.78 c=0.75533091405029995    #corrected3 iteration11
#gamma=0.79 c=0.75116396169710087    #corrected3 iteration11
#gamma=0.8 c=0.74692392576962918    #corrected3 iteration11


theta = 2.0
#alpha = 8.0    # V(x)=rho*x+alpha*x**(1.0/alpha), alpha>1
rho= 2.0    # V(x)=rho*x  
c = 0.46    # from self-consistent calculation using jouwkowsky_parameters_c_selfconsistent_calculation_CR_problem.py for alpha-log(x) potential.
print c

contr = numpy.loadtxt("input_output_files/rho_2.0/contour/nu_contour_c=1_theta=2.0_18000points.txt",float)
# Note that the contour does not depend on parameter 'c' (refer notebook) so we will use the above contour for all different c's.
while True:
        
    def Jc(z):    # joukowsky transformation for hard edge
            return (c*(z+1.0)*(((z+1.0)/z)**(1.0/theta)))
            
    def deriv_Jc(z):    # derivative of joukowsky transformation w.r.t. c
            return (z+1.0)*(((z+1.0)/z)**(1.0/theta))

#    def f1(z):    # See notebook
#            return (1.0/(2.0*(math.pi)*1j))*((b0+b1*Jc(z)+b2*(Jc(z))**2.0+b3*(Jc(z))**3.0+b4*(Jc(z))**4.0)/(z))            

#    def f1(z):    # See notebook
#            return (1.0/(2.0*(math.pi)*1j))*((b0+b1*Jc(z)+b2*(Jc(z))**2+b3*(Jc(z))**3+b4*(Jc(z))**4+b5*(Jc(z))**5+b6*(Jc(z))**6)/(z))

#    def f1(z):    # See notebook
#            return (1.0/(2.0*(math.pi)*1j))*((b0+b1*Jc(z)+b2*(Jc(z))**2+b3*(Jc(z))**3+b4*(Jc(z))**4+b5*(Jc(z))**5+b6*(Jc(z))**6+b7*(Jc(z))**7+b8*(Jc(z))**8)/(z))
            
    def f1(z):    # See notebook
            if (Jc(z)).real<1.0:           
                return (1.0/(2.0*(math.pi)*1j))*((b0+b1*Jc(z)+b2*(Jc(z))**2+b3*(Jc(z))**3+b4*(Jc(z))**4+b5*(Jc(z))**5+b6*(Jc(z))**6+b7*(Jc(z))**7+b8*(Jc(z))**8)/(z))
            else:
                return (1.0/(2.0*(math.pi)*1j))*((c0+c1*Jc(z))/(z))
                
                
    F1 = ((contour_integral.contr_intgl(contr,f1)).real) - (1+theta)
    
#    def deriv_f1(z):    # derivative of f1(z) w.r.t. c
#            return (1.0/(2.0*(math.pi)*1j))*((b1*deriv_Jc(z)+2.0*b2*Jc(z)*deriv_Jc(z)+3.0*b3*((Jc(z))**2.0)*deriv_Jc(z)+4.0*b4*((Jc(z))**3.0)*deriv_Jc(z))/(z))  
                           
    
#    def deriv_f1(z):    # derivative of f1(z) w.r.t. c
#            return (1.0/(2.0*(math.pi)*1j))*((b1*deriv_Jc(z)+2.0*b2*Jc(z)*deriv_Jc(z)+3.0*b3*((Jc(z))**2)*deriv_Jc(z)+4.0*b4*((Jc(z))**3)*deriv_Jc(z)+5.0*b5*((Jc(z))**4)*deriv_Jc(z)+6.0*b6*((Jc(z))**5)*deriv_Jc(z))/(z))    

#    def deriv_f1(z):    # derivative of f1(z) w.r.t. c
#            return (1.0/(2.0*(math.pi)*1j))*((b1*deriv_Jc(z)+2.0*b2*Jc(z)*deriv_Jc(z)+3.0*b3*((Jc(z))**2)*deriv_Jc(z)+4.0*b4*((Jc(z))**3)*deriv_Jc(z)+5.0*b5*((Jc(z))**4)*deriv_Jc(z)+6.0*b6*((Jc(z))**5)*deriv_Jc(z)+7.0*b7*((Jc(z))**6)*deriv_Jc(z)+8.0*b8*((Jc(z))**7)*deriv_Jc(z))/(z))    

    def deriv_f1(z):    # derivative of f1(z) w.r.t. c
            if (Jc(z)).real<1.0:
                return (1.0/(2.0*(math.pi)*1j))*((b1*deriv_Jc(z)+2.0*b2*Jc(z)*deriv_Jc(z)+3.0*b3*((Jc(z))**2)*deriv_Jc(z)+4.0*b4*((Jc(z))**3)*deriv_Jc(z)+5.0*b5*((Jc(z))**4)*deriv_Jc(z)+6.0*b6*((Jc(z))**5)*deriv_Jc(z)+7.0*b7*((Jc(z))**6)*deriv_Jc(z)+8.0*b8*((Jc(z))**7)*deriv_Jc(z))/(z))    
            else:
                return (1.0/(2.0*(math.pi)*1j))*((c1*deriv_Jc(z))/(z))    

    deriv_F1 = ((contour_integral.contr_intgl(contr,deriv_f1)).real)                         
    
    c_next = c - (F1/deriv_F1)    
         
    error = c_next - c
    print('error is', error)
    
    c = c_next

    
#    if(numpy.amax(error)<=1e-2):
    if(abs(error)<=(1e-10)):
        break

print('renormalized parameter c is', c)            