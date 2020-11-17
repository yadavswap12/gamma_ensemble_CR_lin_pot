# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 16:43:33 2018

@author: yadav
"""

#code for computation of parameters c1,c0 of Joukowsky transform for CR paper self consistently. See notebook for details
import math
import cmath
import numpy
import contour_integral
import contour_function
tau = 1.0
rho = -4.0    # V(x)=(tau*x**2 + rho*x)
theta = 2.0
c1 = 69.4   # initializatio for c1(>0), analytical solution is c1=(-2.0/rho), see page number 2439 in C.R.
c0 = 109008.9    # initializatio for c0(>0), analytical solution is c0=(-rho/2.0), see page number 2439 in C.R.        
jacobian = numpy.empty([2,2],float)
ff = numpy.empty(2,float)
c = numpy.empty(2,float)
c[0] = c1
c[1] = c0

while True:
    
    contr = contour_function.contour_CR_prob_soft_edge(c[0],c[1],theta)
#    contr = contour_function.contour_CR_prob_soft_edge(c[0],c[1],theta)    # to make sure output file in contour_function is written properly    
#    contr = contour_function.contour_CR_prob_soft_edge(c[0],c[1],theta)    # to make sure output file in contour_function is written properly    
    
    def f1(z):    # for J(1,1), See notebook
            return (1.0/(2.0*(math.pi)*1j))*((4.0*tau*(c[0]*z+c[1])*z*(((z+1)/z)**(2.0/theta))+rho*z*(((z+1.0)/z)**(1.0/theta)))/(z))
            
    jacobian[0,0] = (contour_integral.contr_intgl(contr,f1)).real 
    print ('J(1,1) is', jacobian[0,0])
       
    def f2(z):    # for J(1,2), See notebook
            return (1.0/(2.0*(math.pi)*1j))*((4.0*tau*(c[0]*z+c[1])*(((z+1.0)/z)**(2.0/theta))+rho*(((z+1.0)/z)**(1.0/theta)))/(z))

    jacobian[0,1] = (contour_integral.contr_intgl(contr,f2)).real 
    print ('J(1,2) is', jacobian[0,1])            
            
    def f3(z):    # for J(2,1), See notebook
            return (1.0/(2.0*(math.pi)*1j))*((4.0*tau*(c[0]*z+c[1])*z*(((z+1.0)/z)**(2.0/theta))+rho*z*(((z+1.0)/z)**(1.0/theta)))/(z+1.0))
            
    jacobian[1,0] = (contour_integral.contr_intgl(contr,f3)).real 
    print ('J(2,1) is', jacobian[1,0])        
            
    def f4(z):    # for J(2,2), See notebook
            return (1.0/(2.0*(math.pi)*1j))*((4.0*tau*(c[0]*z+c[1])*(((z+1.0)/z)**(2.0/theta))+rho*(((z+1.0)/z)**(1.0/theta)))/(z+1.0))
            
    jacobian[1,1] = (contour_integral.contr_intgl(contr,f4)).real 
    print ('J(2,2) is', jacobian[1,1])

    def f5(z):    # for ff1, See notebook
            return (1.0/(2.0*(math.pi)*1j))*((2.0*tau*((c[0]*z+c[1])**2.0)*(((z+1.0)/z)**(2.0/theta))+rho*(c[0]*z+c[1])*(((z+1.0)/z)**(1.0/theta)))/(z)) 

    ff[0] = (contour_integral.contr_intgl(contr,f5)).real  - (1.0+theta) 

    def f6(z):    # for ff2, See notebook
            return (1.0/(2.0*(math.pi)*1j))*((2.0*tau*((c[0]*z+c[1])**2.0)*(((z+1.0)/z)**(2.0/theta))+rho*(c[0]*z+c[1])*(((z+1.0)/z)**(1.0/theta)))/(z+1.0)) 

    ff[1] = (contour_integral.contr_intgl(contr,f6)).real  - 1.0

    delta_X = numpy.linalg.solve(jacobian,ff)            
    
    c_next = c - delta_X    
         
    error = c_next - c
    print('error is', error)
    
    c = c_next
    
    print c[0],c[1]
    
#    if(numpy.amax(error)<=1e-2):
    if(abs(error[0])<=(1e-4) and abs(error[1])<=1e-4 ):
        break            