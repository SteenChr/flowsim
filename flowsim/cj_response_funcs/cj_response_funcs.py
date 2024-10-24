# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 17:16:12 2019

@author: Steen Christensen, Geoscience, AU.
"""
import numpy as np
# import scipy.integrate as integrate
from math import erfc, sqrt, exp, pi, erf, cos, sin, log
from scipy.special import erfcx
from scipy.optimize import brentq
import matplotlib.pyplot as plt
# from time import time

######################################################################################
# Needed math functions
######################################################################################
# Functions to compute roots of x*tan(x) = b and x*cot(x) = b
# Copied from https://stackoverflow.com/questions/53160468/zeros-of-x-tanx-b

def xtan_root(b, n):    # returns the nth root of x*tan(x) - b, starting with n = 0
    if b >= 0:
        left = np.pi*n
        right = left + np.pi/2 - 1e-12
    else:
        right = np.pi*(n + 1)
        left = right - np.pi/2 + 1e-12
    return brentq(lambda x: x*np.tan(x) - b, left, right)

def xcot_root(b, n):    # returns the nth root of x*cot(x) - b, starting with n = 0
    if b <= 1:
        left = np.pi*n + 1e-12
        right = np.pi*(n + 1) - 1e-12
    else:
        left = np.pi*(n + 1) + 1e-12
        right = np.pi*(n + 2) - 1e-12
    return brentq(lambda x: x/np.tan(x) - b, left, right)

######################################################################################
# Response functions from "Conduction of heat in solids" by Carslaw and Jaeger 2nd ed.
######################################################################################
#   App II. Derivatives and integrals of error functions

#   Eq 8. 1st and 2nd devivatives of erf and erfc
def derf(x):
    derf = 2.*exp(-x**2)/sqrt(pi)
    return(derf)
    
def d2erf(x):
    derf2 = -4.*x*exp(-x**2)
    return(derf2)
    
def derfc(x):
    return(-derf(x))
    
def d2erfc(x):
    return(-d2erf(x))
    
#   Eq. 11. First integral
def ierfc(x):
    ierfc = exp(-x*x)/sqrt(pi)
    ierfc -= x*erfc(x)
    return(ierfc)
def dierfc(x):
    dierfc = -2.*x*exp(-x**2)/sqrt(pi) - erfc(x) - x*derfc(x)
    return(dierfc)

#   Eq. 12. Second repeated integral
def i2erfc(x):
    i2erfc = 0.25*(erfc(x) - 2.*x*ierfc(x))
    return(i2erfc)
def di2erfc(x):
    di2erfc = 0.25*(derfc(x) - 2.*(ierfc(x)+x*dierfc(x)))
    return(di2erfc)

def u_erfc(x):
    return(np.frompyfunc(erfc,1,1))

######################################################################################
# Linear reservoir 
######################################################################################
#    
#    Returns a vector of time-dependent response values
#
def lin_res_response(it0,nt,t,x,*args):
    TC =args[0]
    z = np.exp(-np.array(t)/TC)
    z[1:nt] = z[0:nt-1]
    z[0] = 0.0
    return(z)

def lin_res_sens_TC(it0,nt,t,x,*args):
    TC =args[0]
#    z = lin_res_response(it0,nt,t,x,*args)/TC**2
    z = np.exp(-np.array(t)/TC)/TC**2
    z = np.multiply(z,t)
    z[1:nt] = z[0:nt-1]
    z[0] = 0.0

# #   Finite difference approx.
#     fac = 1.00001
#     z1 = lin_res_response(it0,nt,t,x,(TC))
#     z2 = lin_res_response(it0,nt,t,x,(fac*TC))
#     z3 =(z2 - z1)/((fac-1.0)*TC)
#     for it in range(0,20):#nt):
#             print(it,z[it],z3[it],z[it]-z3[it])

    return(z)

    
######################################################################################
# Infinite and semi-infinite solid functions, Ch. 2
######################################################################################
#    Ch. 2.5. Semi-infinite solid functions
#    
#    Returns z which is a vector of time-dependent response values
#
#   Eq. 2: response to b.c. step change
def sinf_solid_step_response(it0,nt,t,x,z0,*args):#,ix,f):
    z = sinf_solid_unit_step_response(it0,nt,t,x,*args)
    return(z0*z)
    
def sinf_solid_unit_step_response(it0,nt,t,x,*args): #D): #,ix,f):
    D = args[0]
    # if isinstance(args, (tuple,list)):
    #     D = args[0]
    # else:
    #     D = arg
    z = np.ones(nt, dtype=np.float64)
    if x != 0.0:
        t0 = t[it0]
        for it in range(it0+1,nt):
            t1 = t[it]-t0
            z[it] = erfc(x/(2.0*sqrt(D*t1)))
#            
#    z1 = np.ones(nt, dtype=np.float64)
#    t0 = t[it0]
#    t1 = t[it0+1:nt] - t0
#    u_erfc = np.frompyfunc(erfc,1,1)
#    y = x/(2.0*np.sqrt(D*t1))
#    z1[it0+1:nt] = u_erfc(y)
            
    #else:
    #    z[it0+1:nt] = 1.0
    return(z)
    
def sinf_solid_unit_step_sensitivity_D(it0,nt,t,x,*args): #D): #,ix,f):
    D = args[0]
    # if isinstance(args, (tuple,list)):
    #     D = args[0]
    # else:
    #     D = args
    z = np.zeros(nt, dtype=np.float64)
    if x != 0.0:
#        xD = x/D
        t0 = t[it0]
        for it in range(it0+1,nt):
            Dt = D*(t[it]-t0)
#            z[it] = 0.5*exp(-0.25*x**2/Dt)*xD/sqrt(pi*Dt)
            b = 0.5*x/sqrt(Dt)
            dbdD = -0.5*b/D
            z[it] = dbdD*derfc(b)
#            print(z[it],z0)
    else:
        z[it0+1:nt] = 0.0
    return(z)
    
def sinf_solid_unit_step_sensitivity_x(it0,nt,t,x,*args): #D): #,ix,f):
    D = args[0]
    # if isinstance(arg, (tuple,list)):
    #     D = args[0]
    # else:
    #     D = args
    z = np.zeros(nt, dtype=np.float64)
#    if x != 0.0:
##        xD = x/D
    t0 = t[it0]
    for it in range(it0+1,nt):
        Dt = D*(t[it]-t0)
#            z[it] = 0.5*exp(-0.25*x**2/Dt)*xD/sqrt(pi*Dt)
        dbdx = 0.5/sqrt(Dt)
        b = dbdx*x
        z[it] = dbdx*derfc(b)
#            print(z[it],z0)
#    else:
#        z[it0+1:nt] = 0.0
    return(z)
    
#   Eqs. 2 and 3: response to temporary b.c. step change (lasts first time interval)
def sinf_solid_temp_step_response(it0,nt,t,x,z0,*args):#,ix,f):
    z = sinf_solid_step_response(it0,nt,t,x,z0,*args)#,ix,f)
    it1 = it0 + 1
    if it1 < nt: 
        z += sinf_solid_step_response(it1,nt,t,x,-z0,*args)#,ix,f)
    return(z)

def sinf_solid_temp_unit_step_response(it0,nt,t,x,*args): #,ix,f):
# This function requires that the time steps are of same length
    z = sinf_solid_unit_step_response(it0,nt,t,x,*args)#,ix,f)
    z1 = np.zeros(nt, dtype=np.float64)
    z1[:] = z[:]
    it1 = it0 + 1
    for it in range(it1+1,nt):
        z[it] = z1[it] - z1[it-1]
    return(z)
    

# Eq. 5: response to linear increase in b.c. value
def sinf_solid_line_response(it0,nt,t,x,k,*args):
    z = sinf_solid_unity_line_response(it0,nt,t,x,*args)
    return(k*z)
#    z = np.zeros(nt, dtype=np.float64)
#    if x != 0.0:
#        t0 = t[it0-1]
#        for it in range(it0,nt):
#            t1 = t[it]-t0
#            z[it] = 4.0*k*t1*i2erfc(x/(2.0*sqrt(D*t1)))
#    return(z)

def sinf_solid_unity_line_response(it0,nt,t,x,*args):
    D = args[0]
    # if isinstance(arg, (tuple,list)):
    #     D = args[0]
    # else:
    #     D = args
    z = np.zeros(nt, dtype=np.float64)
    if x != 0.0:
        t0 = t[it0-1]
        for it in range(it0,nt):
            t1 = t[it]-t0
            z[it] = 4.0*t1*i2erfc(x/(2.0*sqrt(D*t1)))
    else:
        z[it0:nt] = 1.0
    return(z)

#   Eq. 5: response to temporary linear increase in b.c. value (lasts first time interval)
def sinf_solid_temp_linear_response(it0,nt,t,x,k,*args):
    z = sinf_solid_line_response(it0,nt,t,x,k,*args)
    it1 = it0 + 1
    if it1 < nt: 
        z += sinf_solid_line_response(it1,nt,t,x,-k,*args)
    return(z)

def sinf_solid_temp_unity_line_response(it0,nt,t,x,*args):
# This function requires that the time steps are of same length
    z = sinf_solid_unity_line_response(it0,nt,t,x,*args)
    z1 = np.zeros(nt, dtype=np.float64)
    z1[:] = z[:]
    it1 = it0 + 1
    for it in range(it1,nt):
        z[it] = z1[it] - z1[it-1]
    return(z)

######################################################################################
#   Ch. 2.6. Semi-infinite solid harmonic function response
#
# Eq. 8: long-term response harmonic function A*cos(omega*t-epsilon)
def sinf_solid_harmonic_response(nt,t,x,D,A,omega,epsilon):
    z = np.zeros(nt, dtype=np.float64)
    kx=x*sqrt(omega/(2.*D))
    for it in range(0,nt):
        z[it] = A*exp(-kx)*cos(omega*t[it]-kx-epsilon)
    return(z)
    

######################################################################################
#   Ch. 2.8. Semi-infinite solid, radiation into a medium at x = 0
#
# Eq. 2: Radiation into a medium at temperature f(t) = 1 for t > 0
    
def sinf_solid_radiation_unit_response_function(it0,nt,t,x,*args):
    z = np.zeros(nt, dtype=np.float64)
    D = args[0] # Diffusivity
    h = args[1] # Relative boundary conductivity; h = H/K cf. Carslaw and Jaeger, Ch. 1.9 eq. (3)
    t0 = t[it0]
    for it in range(1,nt):
        t1 = t[it]-t0
        a = sqrt(D * t1)
        b = 0.5 * x / a
#        c = h * (x + h * a * a)
        d = b + h*a
        e = erfcx(d) # erfc(x) = exp(-x**2) * erfcx(x), to avoid arithmetic underflow
#        f = c - d**2 + log(e)
        f = log(e) - b**2
        z[it] = erfc(b) - exp(f) # f = exp(c) * erfc(b + h * a)
    # if isinstance(args, (tuple,list)):
    #     D = args[0] # Diffusivity
    #     h = args[1] # Relative boundary conductivity; h = H/K cf. Carslaw and Jaeger, Ch. 1.9 eq. (3)
    #     t0 = t[it0]
    #     for it in range(1,nt):
    #         t1 = t[it]-t0
    #         a = sqrt(D * t1)
    #         b = 0.5 * x / a
    #         c = h * (x + h * a * a)
    #         d = b + h*a
    #         e = erfcx(d) # erfc(x) = exp(-x**2) * erfcx(x), to avoid arithmetic underflow
    #         f = c - d**2 + log(e)
    #         z[it] = erfc(b) - exp(f) # f = exp(c) * erfc(b + h * a)
    # else:
    #     print("*** ERROR *** For solid radiation function, parameter values must be given in tupple!")

    return(z)
    
def sinf_solid_radiation_unit_sensitivity_D(it0,nt,t,x,*args):
    z = np.zeros(nt, dtype=np.float64)
    D = args[0] # Diffusivity
    h = args[1] # Relative boundary conductivity; h = H/K cf. Carslaw and Jaeger, Ch. 1.9 eq. (3)
    t0 = t[it0]
    for it in range(1,nt):#nt):
        t1 = t[it]-t0
        a = sqrt(D * t1)
        b = 0.5 * x / a
        b2 = b**2
#        c = h * (x + h * a * a)
        d = b + h*a
        e = erfcx(d) # erfc(x) = exp(-x**2) * erfcx(x), to avoid arithmetic underflow
        dadD = 0.5*a/D
        dbdD = -0.5*b/D
#            dcdD = t1 * h**2
        dddD = dbdD + h * dadD
#            z[it] = derfc(b) * dbdD - dcdD * exp(c) * erfc(d) - exp(c) * derfc(d) * dddD
#        f1 = c - d**2 + log(e) + 2.*log(h) + log(t1)
        f1 = -b2 + log(e) + 2.*log(h) + log(t1)
#            print(c,-d**2,4.*d/sqrt(pi),dddD)
        sgn = np.sign(dddD)
#        f2 = c - d**2 + log(2./sqrt(pi)) + log(abs(dddD))
        f2 = -b2 + log(2./sqrt(pi)) + log(abs(dddD))
        z[it] = derfc(b) * dbdD - exp(f1) + sgn*exp(f2)# exp(c)*derfc(d)*dddD

##       Finite difference approx.
#        fac = 1.00001
#        z1 = sinf_solid_radiation_unit_response_function(it0,nt,t,x,(D,h))
#        z2 = sinf_solid_radiation_unit_response_function(it0,nt,t,x,(fac*D,h))
#        z3 =(z2 - z1)/((fac-1.0)*D)
#        for it in range(0,nt):
#            print(it,z[it],z3[it],z[it]-z3[it])
#     if isinstance(args, (tuple,list)):
#         D = args[0] # Diffusivity
#         h = args[1] # Relative boundary conductivity; h = H/K cf. Carslaw and Jaeger, Ch. 1.9 eq. (3)
#         t0 = t[it0]
#         for it in range(1,nt):#nt):
#             t1 = t[it]-t0
#             a = sqrt(D * t1)
#             b = 0.5 * x / a
#             c = h * (x + h * a * a)
#             d = b + h*a
#             e = erfcx(d) # erfc(x) = exp(-x**2) * erfcx(x), to avoid arithmetic underflow
#             dadD = 0.5*a/D
#             dbdD = -0.5*b/D
# #            dcdD = t1 * h**2
#             dddD = dbdD + h * dadD
# #            z[it] = derfc(b) * dbdD - dcdD * exp(c) * erfc(d) - exp(c) * derfc(d) * dddD
#             f1 = c - d**2 + log(e) + 2.*log(h) + log(t1)
# #            print(c,-d**2,4.*d/sqrt(pi),dddD)
#             sgn = np.sign(dddD)
#             f2 = c - d**2 + log(2./sqrt(pi)) + log(abs(dddD))
#             z[it] = derfc(b) * dbdD - exp(f1) + sgn*exp(f2)# exp(c)*derfc(d)*dddD

# ##       Finite difference approx.
# #        fac = 1.00001
# #        z1 = sinf_solid_radiation_unit_response_function(it0,nt,t,x,(D,h))
# #        z2 = sinf_solid_radiation_unit_response_function(it0,nt,t,x,(fac*D,h))
# #        z3 =(z2 - z1)/((fac-1.0)*D)
# #        for it in range(0,nt):
# #            print(it,z[it],z3[it],z[it]-z3[it])
#     else:
#         print("*** ERROR *** For solid radiation function, parameter values must be given in tupple!")

    return(z)
    

def sinf_solid_radiation_unit_sensitivity_h(it0,nt,t,x,*args):
    z = np.zeros(nt, dtype=np.float64)
    # if isinstance(args, (tuple,list)):
    D = args[0] # Diffusivity
    h = args[1] # Relative boundary conductivity; h = H/K cf. Carslaw and Jaeger, Ch. 1.9 eq. (3)
    t0 = t[it0]
    for it in range(1,nt):
        t1 = t[it]-t0
        Dt = D * t1
        a = sqrt(Dt)
        b = 0.5 * x / a
        b2 = b**2
#        c = h * (x + h * Dt)
        d = b + h * a
        e = erfcx(d) # erfc(x) = exp(-x**2) * erfcx(x), to avoid arithmetic underflow
#            dcdh = x + 2.* h * Dt
        dddh = a
#        f1 = c - d**2 + log(e) + log(x + 2.*h*Dt)
#        f2 = c - d**2 + log(2./sqrt(pi)) + log(dddh)
        f1 = -b2 + log(e) + log(x + 2.*h*Dt)
        f2 = -b2 + log(2./sqrt(pi)) + log(dddh)
        z[it] = - exp(f1) + exp(f2)# exp(c)*derfc(d)*dddh

##       Finite difference approx.
#        fac = 1.000001
#        z1 = sinf_solid_radiation_unit_response_function(it0,nt,t,x,(D,h))
#        z2 = sinf_solid_radiation_unit_response_function(it0,nt,t,x,(D,fac*h))
#        z3 =(z2 - z1)/((fac-1.0)*h)
#        for it in range(0,nt):
#            print(it,z[it],z3[it],z[it]-z3[it])
        
    # else:
    #     print("*** ERROR *** For solid radiation function, parameter values must be given in tupple!")
    return(z)
    

def sinf_solid_radiation_unit_sensitivity_x(it0,nt,t,x,*args):
    z = np.zeros(nt, dtype=np.float64)
    # if isinstance(args, (tuple,list)):
    D = args[0] # Diffusivity
    h = args[1] # Relative boundary conductivity; h = H/K cf. Carslaw and Jaeger, Ch. 1.9 eq. (3)
    t0 = t[it0]
    for it in range(1,nt):#nt):
        t1 = t[it]-t0
        a = sqrt(D * t1)
        dbdx = 0.5 / a
        b = dbdx * x
        b2 = b**2
#        c = h * (x + h * a * a)
#            dcdx = c - h*x + h
        d = b + h*a
        dddx = dbdx
        e = erfcx(d) # erfc(x) = exp(-x**2) * erfcx(x), to avoid arithmetic underflow
#            z[it] = derfc(b) * dbdx - dcdx * exp(c) * erfc(d) - exp(c) * derfc(d) * dddx
#        f1 = c - d**2 + log(e) + log(h) # log(dcdx*exp(c)*erfc(d)) = log(dcdx*exp(c)*exp(-d**2)*erfcx(d))
#        f2 = c - d**2 + log(2./sqrt(pi)) + log(dddx) # log(-exp(c)*derfc(d)*dddx)
        f1 = -b2 + log(e) + log(h) # log(dcdx*exp(c)*erfc(d)) = log(dcdx*exp(c)*exp(-d**2)*erfcx(d))
        f2 = -b2 + log(2./sqrt(pi)) + log(dddx) # log(-exp(c)*derfc(d)*dddx)
        z[it] = derfc(b) * dbdx - exp(f1) + exp(f2)
#
##       Finite difference approx.
#        fac = 1.00001
#        z1 = sinf_solid_radiation_unit_response_function(it0,nt,t,x,(D,h))
#        z2 = sinf_solid_radiation_unit_response_function(it0,nt,t,fac*x,(D,h))
#        z3 =(z2 - z1)/((fac-1.0)*x)
#        for it in range(0,nt):
#            print(it,z[it],z3[it],z[it]-z3[it])
    # else:
    #     print("*** ERROR *** For solid radiation function, parameter values must be given in tupple!")

    return(z)
   

######################################################################################
#   Ch. 2.11. Semi-infinite solid with heat generated within it
#
# Eq. 2: Heat production at constant rate (A = 1), with a = b = 0 (initial temp. = 0)
    
def sinf_solid_unit_production_response_function(it0,nt,t,x,*args):
    z = np.zeros(nt, dtype=np.float64)
    # if isinstance(args, (tuple,list)):
    D = args[0] # Diffusivity
    K = args[1] # Conductivity
    if x != 0.0:
        t0 = t[it0]
        c1 = 0.5*(x**2)/K
        for it in range(it0+1,nt):
            Dt = D*(t[it]-t0)
            c2 = 0.5*x/sqrt(Dt)
            c3 = sqrt(Dt/pi)*x/K
            z[it] = (Dt/K + c1)*erf(c2) + c3*exp(-(c2**2)) - c1
    # else:
    #     print("*** ERROR *** For solid production function, parameter values must be given in tupple!")
    return(z)

def sinf_solid_unit_production_sensitivity_D(it0,nt,t,x,*args):
    z = np.zeros(nt, dtype=np.float64)
    # if isinstance(args, (tuple,list)):
    D = args[0] # Diffusivity
    K = args[1] # Conductivity
    if x != 0.0:
        t0 = t[it0]
        for it in range(it0+1,nt):
            delt = t[it]-t0
            Dt = D*delt
            c2 = 0.5*x/sqrt(Dt)
            z[it] = erf(c2)*delt/K
#                
#                a = (Dt + 0.5*x**2)/K
#                b = 0.5*x/sqrt(Dt)
#                c = x*sqrt(Dt/pi)/K
#                d = -b**2
#                dadD = delt/K
#                dbdD = -0.5*b/D
#                dcdD = 0.5*c/D
#                dddD = -2.*b*dbdD
#                z0  = dadD * erf(b) + a * derf(b) * dbdD + (dcdD + c * dddD) * exp(d)
#
#                z1 = a*erf(b)+c*exp(d)-0.5*x**2/K
#                fac=1.0001
#                delD=(fac-1.)*D
#                D=fac*D
#                Dt = D*delt
#                a = (Dt + 0.5*x**2)/K
#                b = 0.5*x/sqrt(Dt)
#                c = x*sqrt(Dt/pi)/K
#                d = -b**2
#                z2 = a*erf(b)+c*exp(d)-0.5*x**2/K
#                print(z[it], z0, (z2-z1)/delD)
#                print(dadD * erf(b), a*derf(b)*dbdD, dcdD*exp(d), c*dddD*exp(d))
    # else:
    #     print("*** ERROR *** For solid production sensitivity function, parameter values must be given in tupple!")
    return(z)

def sinf_solid_unit_production_sensitivity_K(it0,nt,t,x,*args):
    z = np.zeros(nt, dtype=np.float64)
    # if isinstance(args, (tuple,list)):
#    D = args[0] # Diffusivity
    K = args[1] # Conductivity
    z = sinf_solid_unit_production_response_function(it0,nt,t,x,*args)        
    # else:
    #     print("*** ERROR *** For solid production sensitivity function, parameter values must be given in tupple!")
    return(-z/K)
    
def sinf_solid_unit_production_sensitivity_x(it0,nt,t,x,*args):
    z = np.zeros(nt, dtype=np.float64)
    # if isinstance(args, (tuple,list)):
    D = args[0] # Diffusivity
    K = args[1] # Conductivity
#        if x != 0.0:
    t0 = t[it0]
    c1d = x/K
    c1 = 0.5*x*c1d
    for it in range(it0+1,nt):
        Dt = D*(t[it]-t0)
        c2d = 0.5/sqrt(Dt)
        c2 = c2d*x
        c3d = sqrt(Dt/pi)/K
        c3 = c3d*x
        z[it] = ((Dt/K + c1)*c2d*derf(c2) +
                 (c3d - c3*2.*c2*c2d)*exp(-(c2**2)) -
                 c1d + c1d*erf(c2))
    # else:
    #     print("*** ERROR *** For solid production sensitivity function, parameter values must be given in tupple!")

# #   Finite difference approx.
#     fac = 1.00001
#     z1 = sinf_solid_unit_production_response_function(it0,nt,t,x,*(D,K))
#     z2 = sinf_solid_unit_production_response_function(it0,nt,t,x+.00001,*(D,K))
#     z3 =(z2 - z1)/.00001#((fac-1.0)*x)
#     for it in range(0,nt):
#         print(it,z[it],z3[it],z[it]-z3[it])

    return(z)

######################################################################################
"""
Ch. 12.4. Semi-infinite solid

Eq. 23: Heat production at constant rate (A = 1; or k=1 and n=0),
radiation into medium at zero temperature, solids initial temp. = 0
"""
    
def sinf_solid_radiation_unit_production_response_function(it0,nt,t,x,*args):
    z = np.zeros(nt, dtype=np.float64)
    D = args[0] # Diffusivity
    K = args[1] # Conductivity
    h = args[2]  # Relative boundary conductivity; h = H/K cf. Carslaw and Jaeger, Ch. 1.9 eq. (3)
    t0 = t[it0]
    for it in range(it0+1,nt):
        a = sqrt(D*(t[it]-t0))
        b = 0.5*x/a
        d = b + h*a
        e = erfcx(d)
        f = exp(log(e) - b**2)
        g = erfc(b) - 2.*h*a*ierfc(b) + 4.*i2erfc(b)*(h*a)**2
        z[it] =a**2 + (f - g)/h**2
    z = z/K
    return(z)

def sinf_solid_radiation_unit_production_sensitivity_D(it0,nt,t,x,*args):
    z = np.zeros(nt, dtype=np.float64)
    D = args[0] # Diffusivity
    K = args[1] # Conductivity
    h = args[2]  # Relative boundary conductivity; h = H/K cf. Carslaw and Jaeger, Ch. 1.9 eq. (3)
    t0 = t[it0]
    for it in range(it0+1,nt):
        t1 = t[it]-t0
        a = sqrt(D*t1)
        dadD = 0.5*sqrt(t1/D)
        b = 0.5*x/a
        dbdD = -0.5*b/D
        dcdD = t1*h**2
        d = b + h*a
        dddD = dbdD + h*dadD
        e = erfcx(d)
        f = dcdD*exp(log(e) - b**2) - 2.*dddD*exp(-b**2)/sqrt(pi)
        g = dbdD*(derfc(b) - 2.*h*a*dierfc(b) + 4.*di2erfc(b)*(h*a)**2)
        g = g - 2.*h*dadD*ierfc(b) + 8.*(h**2)*a*dadD*i2erfc(b)
        z[it] = (f - g)/(K*h**2) + t1/K
# #       Finite difference approx.
#     fac = 1.000001
#     z1 = sinf_solid_radiation_unit_production_response_function(it0,nt,t,x,*(D,K,h))
#     z2 = sinf_solid_radiation_unit_production_response_function(it0,nt,t,x,*(fac*D,K,h))
#     z3 =(z2 - z1)/((fac-1.0)*D)
#     for it in range(0,nt):
#         print(it,z[it],z3[it],z[it]-z3[it])
    return(z)

def sinf_solid_radiation_unit_production_sensitivity_h(it0,nt,t,x,*args):
    z = np.zeros(nt, dtype=np.float64)
    D = args[0] # Diffusivity
    K = args[1] # Conductivity
    h = args[2]  # Relative boundary conductivity; h = H/K cf. Carslaw and Jaeger, Ch. 1.9 eq. (3)
    t0 = t[it0]
    for it in range(it0+1,nt):
        t1 = t[it]-t0
        a = sqrt(D*t1)
        b = 0.5*x/a
        # c = h*x + (h*a)**2
        d = b + h*a
        e = erfcx(d)
        # dadh = 0.0
        # dbdh = 0.0
        dcdh = x + 2.*h*a**2
        dddh = a
        f = exp(log(e) - b**2)
        f = f - (erfc(b) - 2*h*a*ierfc(b) + 
                  4.*i2erfc(b)*(h*a)**2)
        f = -2.*f/h**3
        g = exp(log(dcdh) + log(e) - b**2) - 2.*exp(log(dddh) - b**2 - log(sqrt(pi)))
        g = g + 2.*a*ierfc(b) - 8.*h*i2erfc(b)*a**2
        g = g/h**2
        z[it] = (f + g)/K
        
# #       Finite difference approx.
#     fac = 1.001
#     z1 = sinf_solid_radiation_unit_production_response_function(it0,nt,t,x,*(D,K,h))
#     z2 = sinf_solid_radiation_unit_production_response_function(it0,nt,t,x,*(D,K,fac*h))
#     z3 =(z2 - z1)/((fac-1.0)*h)
#     for it in range(0,200):#nt):
#         if z3[it]!=0.0:
#             print(it,z[it],z3[it],z[it]-z3[it],(z[it]-z3[it])/z3[it])
    return(z)

def sinf_solid_radiation_unit_production_sensitivity_K(it0,nt,t,x,*args):
    # z = np.zeros(nt, dtype=np.float64)
    D = args[0] # Diffusivity
    K = args[1] # Conductivity
    h = args[2]  # Relative boundary conductivity; h = H/K cf. Carslaw and Jaeger, Ch. 1.9 eq. (3)
    z = -sinf_solid_radiation_unit_production_response_function(it0,nt,t,x,
                                                                *args)/K
    z = z + sinf_solid_radiation_unit_production_sensitivity_D(it0,nt,t,x,
                                                               *args)*D/K
    z = z - sinf_solid_radiation_unit_production_sensitivity_h(it0,nt,t,x,
                                                               *args)*h/K
# #       Finite difference approx.
#     fac = 1.001
#     z1 = sinf_solid_radiation_unit_production_response_function(it0,nt,t,x,*(D,K,h))
#     z2 = sinf_solid_radiation_unit_production_response_function(it0,nt,t,x,
#                                                                 *(fac*D,fac*K,h/fac))
#     z3 =(z2 - z1)/((fac-1.0)*K)
#     for it in range(0,200):#nt):
#         if z3[it]!=0.0:
#             print(it,z[it],z3[it],z[it]-z3[it],(z[it]-z3[it])/z3[it])
    return(z)


def sinf_solid_radiation_unit_production_sensitivity_x(it0,nt,t,x,*args):
    z = np.zeros(nt, dtype=np.float64)
    D = args[0] # Diffusivity
    K = args[1] # Conductivity
    h = args[2]  # Relative boundary conductivity; h = H/K cf. Carslaw and Jaeger, Ch. 1.9 eq. (3)
    t0 = t[it0]
    for it in range(it0+1,nt):
        a = sqrt(D*(t[it]-t0))
        dbdx = 0.5/a
        b = x*dbdx
        d = b + h*a
        e = erfcx(d)
        f = h*exp(log(e) - b**2) - exp(-b**2)/a/sqrt(pi)
        g = dbdx*(derfc(b) - 2.*h*a*dierfc(b) + 4.*di2erfc(b)*(h*a)**2)
        z[it] = f - g
    z = z/(K*h**2)
# #   Finite difference approx.
#     fac = 1.00001
#     z1 = sinf_solid_radiation_unit_production_response_function(it0,nt,t,x,*args)
#     z2 = sinf_solid_radiation_unit_production_response_function(it0,nt,t,x+.001,*args)
#     z3 =(z2 - z1)/.001#((fac-1.0)*x)
#     for it in range(0,200):#nt):
#         print(it,z[it],z3[it],z[it]-z3[it])    
    return(z)

######################################################################################
# Solid slab functions, Ch. 3.
######################################################################################
"""
Ch. 3.12, eq. 3. Slab with zero initial temperature and radiation at ends 
into a medium at temperature f(t) = 1 for t > 0
"""

def slab_radiation_unit_response_function(it0,nt,t,x,*args):
    z = np.zeros(nt, dtype=np.float64)
    # if isinstance(args, (tuple,list)):
    D = args[0] # Diffusivity
    h = args[1] # Relative boundary conductivity; h = H/K cf. Carslaw and Jaeger, Ch. 1.9 eq. (3)
    l = args[2] # Length of slab

#   Use coordinate transform, so z-b.c. wil be at x = 0, not at x = -l, 
#   and zero gradient will be at x = l, not at x = 0.
    xl = x - l

    nmax = 1000
    hl = h*l
    alpha = np.zeros(nmax, dtype=float)
    for n in range(0,nmax):
        alpha[n] = xtan_root(hl, n)/l
    alpha2 = np.power(alpha, 2)
        
    t0 = t[it0]
    h2 = h*h
    EPS = 1.e-12
    for it in range(it0+1,nt):
        at = D*(t[it]-t0)
        z[it] = 1.
        dzsum = 0.0
        dz = 1.
        n = 0
        while np.abs(dz) > EPS and n < nmax:
#            while n < nmax:
            dz = (np.exp(-at*alpha2[n])
                 * (2.*h/((h2+alpha2[n])*l+h))
                 * (np.cos(alpha[n]*xl)/np.cos(alpha[n]*l))
                 )
            dzsum = dzsum + dz
            n += 1
        z[it] = 1. - dzsum
    # else:
    #     print("*** ERROR *** For slab radiation function, parameter values must be given in tupple!")
    return(z)
    

    
def slab_radiation_sensitivity_D(it0,nt,t,x,*args):
    z = np.zeros(nt, dtype=np.float64)
    # if isinstance(args, (tuple,list)):
    D = args[0] # Diffusivity
    h = args[1] # Relative boundary conductivity; h = H/K cf. Carslaw and Jaeger, Ch. 1.9 eq. (3)
    l = args[2] # Length of slab

#   Use coordinate transform, so z-b.c. wil be at x = 0, not at x = -l, 
#   and zero gradient will be at x = l, not at x = 0.
    xl = x - l

    nmax = 1000
    hl = h*l
    alpha = np.zeros(nmax, dtype=float)
    for n in range(0,nmax):
        alpha[n] = xtan_root(hl, n)/l
    alpha2 = np.power(alpha, 2)
        
    t0 = t[it0]
    h2 = h*h
    EPS = 1.e-12
    for it in range(it0+1,nt):
        at = D*(t[it]-t0)
        t1 = t[it] - t0
        z[it] = 1.
        dzsum = 0.0
        dz = 1.
        n = 0
        while np.abs(dz) > EPS and n < nmax:
#            while n < nmax:
            dz = alpha2[n]*t1*(np.exp(-at*alpha2[n])
                 * (2.*h/((h2+alpha2[n])*l+h))
                 * (np.cos(alpha[n]*xl)/np.cos(alpha[n]*l))
                 )
            dzsum = dzsum + dz
            n += 1
        z[it] = dzsum
    # else:
    #     print("*** ERROR *** For slab radiation sensitivity_D function, parameter values must be given in tupple!")
    return(z)
    
def slab_radiation_sensitivity_h(it0,nt,t,x,*args):
    """
    Analytical solution is not derived because of complications.
    Reason: alpha depends on h.
    Finite difference approximation is used instead.
    """
    z = np.zeros(nt, dtype=np.float64)
    D = args[0] # Diffusivity
    h = args[1] # Relative boundary conductivity; h = H/K cf. Carslaw and Jaeger, Ch. 1.9 eq. (3)
    l = args[2] # Length of slab
#   Finite difference approx.
    fac = 1.00001
    z = (slab_radiation_unit_response_function(it0,nt,t,x,*(D,h*fac,l)) -
         slab_radiation_unit_response_function(it0,nt,t,x,*(D,h,l)))
    z = z/((fac-1.0)*h)
    return(z)    
    
def slab_radiation_sensitivity_l(it0,nt,t,x,*args):
    """
    Analytical solution is not derived because of complications.
    Reason: alpha depends on h.
    Finite difference approximation is used instead.
    """
    z = np.zeros(nt, dtype=np.float64)
    D = args[0] # Diffusivity
    h = args[1] # Relative boundary conductivity; h = H/K cf. Carslaw and Jaeger, Ch. 1.9 eq. (3)
    l = args[2] # Length of slab
#   Finite difference approx.
    fac = 1.00001
    z = (slab_radiation_unit_response_function(it0,nt,t,x,*(D,h,l*fac)) -
         slab_radiation_unit_response_function(it0,nt,t,x,*(D,h,l)))
    z = z/((fac-1.0)*l)
    return(z)    
    
def slab_radiation_sensitivity_x(it0,nt,t,x,*args):
    z = np.zeros(nt, dtype=np.float64)
    # if isinstance(args, (tuple,list)):
    D = args[0] # Diffusivity
    h = args[1] # Relative boundary conductivity; h = H/K cf. Carslaw and Jaeger, Ch. 1.9 eq. (3)
    l = args[2] # Length of slab

#   Use coordinate transform, so z-b.c. wil be at x = 0, not at x = -l, 
#  and zero gradient will be at x = l, not at x = 0.
    xl = x - l

    nmax = 1000
    hl = h*l
    alpha = np.zeros(nmax, dtype=float)
    for n in range(0,nmax):
        alpha[n] = xtan_root(hl, n)/l
    alpha2 = np.power(alpha, 2)
        
    t0 = t[it0]
    h2 = h*h
    EPS = 1.e-12
    for it in range(it0+1,nt):
        at = D*(t[it]-t0)
        z[it] = 1.
        dzsum = 0.0
        dz = 1.
        n = 0
        while np.abs(dz) > EPS and n < nmax:
            dz = (np.exp(-at*alpha2[n])
                 * (2.*h/((h2+alpha2[n])*l+h))
                 * (alpha[n]*np.sin(alpha[n]*xl)/np.cos(alpha[n]*l))
                 )
            dzsum = dzsum + dz
            n += 1
        z[it] = dzsum
    # else:
    #     print("*** ERROR *** For slab radiation sensitivity_h function, parameter values must be given in tupple!")
    return(z)

######################################################################################
"""
Ch. 3.5, eq. 6. Slab with with zero initial temperature and 
unit temperature at ends, f(t) = 1, for t > 0.
"""
def sup_term00(n,pixl,at):#,sens_D):
    m = 2*n+1
    a = (-1)**n/m
    b = cos(m*pixl)
    d = -at*m**2
    c = exp(d)
    return(a*b*c)
def slab_unit_step_response(it0,nt,t,x,*args):
    z = np.zeros(nt, dtype=np.float64)
    z[it0+1:nt] = 1.0
    if x > 0.0:
        D = args[0] # Diffusivity
        l = args[1] # Length of slab
#       Use coordinate transform, so z-b.c. wil be at x = 0, not at x = l, 
#       and zero gradient will be at x = l, not at x = 0.
        nmax = 500
        tol = 1.e-8
        xl = x/l - 1.
        pixl = 0.5 * pi * xl
        a = 0.25 * D * (pi/l)**2
        b = 4./pi
        t0 = t[it0]
        for it in range(it0+1,nt):
            at = a*(t[it]-t0)
            dzsum0 = 0.0
            dzsum = 0.0
            for n in range(0,nmax):
                dz = sup_term00(n,pixl,at)#,sens_D)
                z[it] -= b * dz
                dzsum += dz
                if abs(dzsum-dzsum0)<tol:
                    break
                else:
                    dzsum0 = dzsum
    return(z)

def sup_term01(n,pixl,at):#,sens_D):
    m = 2*n+1
    a = m*(-1)**n
    b = cos(m*pixl)
    d = -at*m**2
    c = exp(d)
    return(a*b*c)
def slab_unit_step_sensitivity_D(it0,nt,t,x,*args):
    z = np.zeros(nt, dtype=np.float64)
    if x > 0.0:
        D = args[0] # Diffusivity
        l = args[1] # Length of slab
#       Use coordinate transform, so z-b.c. wil be at x = 0, not at x = l, 
#       and zero gradient will be at x = l, not at x = 0.
        nmax = 500
        tol = 1.e-8
        xl = x/l - 1.
        pixl = 0.5 * pi * xl
        a = 0.25 * D * (pi/l)**2
        b = pi/l**2
        t0 = t[it0]
        for it in range(it0+1,nt):
            t1 = t[it]-t0
            at = a*t1
            dzsum0 = 0.0
            dzsum = 0.0
            for n in range(0,nmax):
                dz = sup_term01(n,pixl,at)#,sens_D)
                z[it] += dz
                dzsum += dz
                if abs(dzsum-dzsum0)<tol:
                    break
                else:
                    dzsum0 = dzsum
            z[it] *= b*t1
        
# #       Finite difference approx.
#     fac = 1.000001
#     z1 = slab_unit_step_response(it0,nt,t,x,*(D,l))
#     z2 = slab_unit_step_response(it0,nt,t,x,*(fac*D,l))
#     z3 =(z2 - z1)/((fac-1.0)*D)
#     for it in range(0,nt):
#         print(it,z[it],z3[it],z[it]-z3[it])
        
    return(z)

def sup_term02(n,pixl,at):#,sens_D):
    m = 2*n+1
    a = (-1)**n
    b = sin(m*pixl)
    d = -at*m**2
    c = exp(d)
    return(a*b*c)
def slab_unit_step_sensitivity_x(it0,nt,t,x,*args):
    z = np.zeros(nt, dtype=np.float64)
    D = args[0] # Diffusivity
    l = args[1] # Length of slab
#   Use coordinate transform, so z-b.c. wil be at x = 0, not at x = l, 
#   and zero gradient will be at x = l, not at x = 0.
#    if x != 0:
    nmax = 500
    tol = 1.e-8
    xl = x/l - 1.
    pixl = 0.5 * pi * xl
    a = 0.25 * D * (pi/l)**2
    b = 2./l
    t0 = t[it0]
    for it in range(it0+1,nt):
        at = a*(t[it]-t0)
        dzsum0 = 0.0
        dzsum = 0.0
        for n in range(0,nmax):
            dz = sup_term02(n,pixl,at)
            z[it] += dz
            dzsum += dz
            if abs(dzsum-dzsum0)<tol:
                break
            else:
                dzsum0 = dzsum
        z[it] *= b
# #   Finite difference approx.
#     fac = 1.00001
#     z1 = slab_unit_step_response(it0,nt,t,x,*args)
#     z2 = slab_unit_step_response(it0,nt,t,x+.00001,*args)
#     z3 =(z2 - z1)/.00001#((fac-1.0)*x)
#     for it in range(0,nt):
#         print(it,z[it],z3[it],z[it]-z3[it])

    return(z)
  
######################################################################################
"""
Ch. 3.14, eq. 7. Slab with heat generated within it. Heat production at 
constant rate (A = 1), and zero temperature at ends.
"""
    
def sup_term(n,pixl,at):#,sens_D):
    m = 2*n+1
    a = (-1)**n/m**3
    b = cos(m*pixl)
    d = -at*m**2
    c = exp(d)
#    if not sens_D:
#        d = 1.0
#    return(a*b*c*d)
    return(a*b*c)
    
def slab_unit_production_response_function(it0,nt,t,x,*args):
    z = np.zeros(nt, dtype=np.float64)
    # if isinstance(args, (tuple,list)):
    D = args[0] # Diffusivity
    K = args[1] # Conductivity
    l = args[2] # Length of slab
#   Use coordinate transform, so z-b.c. wil be at x = 0, not at x = l, 
#   and zero gradient will be at x = l, not at x = 0.
#    if x != 0:
    nmax = 500
    tol = 1.e-8
#        xl = 1. - x/l
    xl = x/l - 1.
    z[it0+1:nt] = 1.0 - xl**2
    pixl = 0.5 * pi * xl
    a = 0.25 * D * (pi/l)**2
    b = 32./pi**3
    t0 = t[it0]
    for it in range(it0+1,nt):
        at = a*(t[it]-t0)
        dzsum0 = 0.0
        dzsum = 0.0
        for n in range(0,nmax):
            dz = sup_term(n,pixl,at)#,sens_D)
            z[it] -= b * dz
            dzsum += dz
            if abs(dzsum-dzsum0)<tol:
                break
            else:
                dzsum0 = dzsum
    c = 0.5 * l**2 / K
    z = c * z
    # else:
    #     print("*** ERROR *** For slab production function, parameter values must be given in tupple!")
    return(z)

#def slab_unit_production_response_function(it0,nt,t,x,*args):
#    z = np.zeros(nt, dtype=np.float64)
#    if isinstance(args, (tuple,list)):
#        sens_D = False
#        z = slab_unit_production_response(it0,nt,t,x,sens_D,*args)
#    else:
#        print("*** ERROR *** For slab production function, parameter values must be given in tupple!")
#    return(z)

def sup_term1(n,pixl,at):#,sens_D):
    m = 2*n+1
    a = (-1)**n/m**3
    b = cos(m*pixl)
    d = -at*m**2
    c = exp(d)
#    if not sens_D:
#        d = 1.0
    return(a*b*c*d)
def slab_unit_production_sensitivity_D(it0,nt,t,x,*args):
    z = np.zeros(nt, dtype=np.float64)
    # if isinstance(args, (tuple,list)):
    D = args[0] # Diffusivity
    K = args[1] # Conductivity
    l = args[2] # Length of slab
#   Use coordinate transform, so z-b.c. wil be at x = 0, not at x = l, 
#   and zero gradient will be at x = l, not at x = 0.
#    if x != 0:
    nmax = 500
    tol = 1.e-8
#    xl = 1. - x/l
    xl = x/l - 1.
#    z[it0+1:nt] = 1.0 - xl**2
    pixl = 0.5 * pi * xl
    a = 0.25 * D * (pi/l)**2
    b = 32./pi**3
    t0 = t[it0]
    for it in range(it0+1,nt):
        at = a*(t[it]-t0)
        dzsum0 = 0.0
        dzsum = 0.0
        for n in range(0,nmax):
            dz = sup_term1(n,pixl,at)#/D#,sens_D)
            z[it] -= dz
            dzsum += dz
            if abs(dzsum-dzsum0)<tol:
                break
            else:
                dzsum0 = dzsum
    c = 0.5 * l**2 *b / K
    z = c * z/D
    # else:
        # print("*** ERROR *** For slab production function, parameter values must be given in tupple!")
    return(z)

#    z = np.zeros(nt, dtype=np.float64)
#    if isinstance(args, (tuple,list)):
#        D = args[0] # Diffusivity
#        sens_D = True
#        z = slab_unit_production_response(it0,nt,t,x,sens_D,*args)
#        z = z/D
#        
##       Finite difference approx.
#        D = args[0] # Diffusivity
#        K = args[1] # Conductivity
#        l = args[2] # Length of slab
#        fac = 1.000001
#        z1 = slab_unit_production_response(it0,nt,t,x,False,*(D,K,l))
#        z2 = slab_unit_production_response(it0,nt,t,x,False,*(fac*D,K,l))
#        z3 =(z2 - z1)/((fac-1.0)*D)
#        for it in range(0,nt):
#            print(it,z[it],z3[it],z[it]-z3[it])
#    else:
#        print("*** ERROR *** For slab production sensitivity function, parameter values must be given in tupple!")
#    return(z)

def slab_unit_production_sensitivity_K(it0,nt,t,x,*args):
    D = args[0] # Diffusivity
    K = args[1] # Conductivity
#    l = args[2] # Length of slab
    z = -slab_unit_production_response_function(it0,nt,t,x,*args)
    z = z + slab_unit_production_sensitivity_D(it0,nt,t,x,*args)*D
    z = z/K
# #       Finite difference approx.
#     fac = 1.000001
#     z1 = slab_unit_production_response_function(it0,nt,t,x,*(D,K,l))
#     z2 = slab_unit_production_response_function(it0,nt,t,x,*(fac*D,fac*K,l))
#     z3 =(z2 - z1)/((fac-1.0)*K)
#     for it in range(0,100):#nt):
#         if z[it] != 0.0:
#             print(it,z[it],z3[it],(z[it]-z3[it])/z3[it])
    return(z)
    
def sup_term2(n,pixl,at):#,sens_D):
    m = 2*n+1
    a = (-1)**n/m**2
    b = sin(m*pixl)
    d = -at*m**2
    c = exp(d)
#    if not sens_D:
#        d = 1.0
#    return(a*b*c*d)
    return(a*b*c)
    
def slab_unit_production_sensitivity_x(it0,nt,t,x,*args):
    z = np.zeros(nt, dtype=np.float64)
    # if isinstance(args, (tuple,list)):
    D = args[0] # Diffusivity
    K = args[1] # Conductivity
    l = args[2] # Length of slab
#   Use coordinate transform, so z-b.c. wil be at x = 0, not at x = l, 
#   and zero gradient will be at x = l, not at x = 0.
#    if x != 0:
    nmax = 500
    tol = 1.e-8
#        xl = 1. - x/l
    xl = x/l - 1.
    z[it0+1:nt] = -2.*xl/l
    pixl = 0.5 * pi * xl
    a = 0.25 * D * (pi/l)**2
    b = 16./pi**2/l
    t0 = t[it0]
    for it in range(it0+1,nt):
        at = a*(t[it]-t0)
        dzsum0 = 0.0
        dzsum = 0.0
        for n in range(0,nmax):
            dz = sup_term2(n,pixl,at)#,sens_D)
            z[it] += b * dz
            dzsum += dz
            if abs(dzsum-dzsum0)<tol:
                break
            else:
                dzsum0 = dzsum
    c = 0.5 * l**2 / K
    z = c * z
    # else:
        # print("*** ERROR *** For slab production function, parameter values must be given in tupple!")
#
##   Finite difference approx.
#    fac = 1.00001
#    z1 = slab_unit_production_response_function(it0,nt,t,x,*args)
#    z2 = slab_unit_production_response_function(it0,nt,t,x+.001,*args)
#    z3 =(z2 - z1)/.001#((fac-1.0)*x)
#    for it in range(0,nt):
#        print(it,z[it],z3[it],z[it]-z3[it])

    return(z)
    
######################################################################################
"""
Ch. 3.14, eq. 12. Slab with zero initial temperature, radiation at ends into
medium at zero temperature, and constant unit production for t>0.
"""

def slab_radiation_unit_production_response_function(it0,nt,t,x,*args):
    z = np.zeros(nt, dtype=np.float64)
    # if isinstance(args, (tuple,list)):
    D = args[0] # Diffusivity
    K = args[1] # Conductivity
    h = args[2] # Relative boundary conductivity; h = H/K cf. Carslaw and Jaeger, Ch. 1.9 eq. (3)
    l = args[3] # Length of slab
    
#   Use coordinate transform, so z-b.c. wil be at x = 0, not at x = -l, 
#   and zero gradient will be at x = l, not at x = 0.
    xl = x - l

    nmax = 1000
    hl = h*l
    alpha = np.zeros(nmax, dtype=float)
    for n in range(0,nmax):
        alpha[n] = xtan_root(hl, n)/l
        if alpha[n] < 0.0:
            print("\nalpha[%d]=%f < 0.0"%(n,alpha[n]))
    alpha2 = np.power(alpha, 2)
        
    t0 = t[it0]
    h2 = h*h
    Kh = 0.5/(K*h)
    h0 = 2.*l + h*(l**2-xl**2)
    EPS = 1.e-12
    for it in range(it0+1,nt):
        at = D*(t[it]-t0)
        z[it] = 1.
        dzsum = 0.0
        dz = 1.
        n = 0
        while np.abs(dz) > EPS and n < nmax:
#            while n < nmax:
            dz = (np.exp(-at*alpha2[n])
                 * (1./(alpha2[n]*((h2+alpha2[n])*l+h)))
                 * (np.cos(alpha[n]*xl)/np.cos(alpha[n]*l))
                 )
            dzsum = dzsum + dz
            n += 1
        z[it] = Kh * (h0 - 4.*h2*dzsum)
    # else:
    #     print("*** ERROR *** For slab radiation/production function, parameter values must be given in tupple!")
    return(z)    

def slab_radiation_unit_production_sensitivity_D(it0,nt,t,x,*args):
    z = np.zeros(nt, dtype=np.float64)
    # if isinstance(args, (tuple,list)):
    D = args[0] # Diffusivity
    K = args[1] # Conductivity
    h = args[2] # Relative boundary conductivity; h = H/K cf. Carslaw and Jaeger, Ch. 1.9 eq. (3)
    l = args[3] # Length of slab
    
#   Use coordinate transform, so z-b.c. wil be at x = 0, not at x = -l, 
#   and zero gradient will be at x = l, not at x = 0.
    xl = x - l

    nmax = 1000
    hl = h*l
    alpha = np.zeros(nmax, dtype=float)
    for n in range(0,nmax):
        alpha[n] = xtan_root(hl, n)/l
        if alpha[n] < 0.0:
            print("\nalpha[%d]=%f < 0.0"%(n,alpha[n]))
    alpha2 = np.power(alpha, 2)
        
    t0 = t[it0]
    h2 = h*h
    Kh = 0.5/(K*h)
    EPS = 1.e-12
    for it in range(it0+1,nt):
        t1 = t[it] - t0
        at = D*t1
        z[it] = 1.
        dzsum = 0.0
        dz = 1.
        n = 0
        while np.abs(dz) > EPS and n < nmax:
#            while n < nmax:
            dz = alpha2[n]*t1*(np.exp(-at*alpha2[n])
                 * (1./(alpha2[n]*((h2+alpha2[n])*l+h)))
                 * (np.cos(alpha[n]*xl)/np.cos(alpha[n]*l))
                 )
            dzsum = dzsum + dz
            n += 1
        z[it] = Kh * (4.*h2*dzsum)
    # else:
    #     print("*** ERROR *** For slab radiation/production sensitivity_D function, parameter values must be given in tupple!")
    return(z)    

def slab_radiation_unit_production_sensitivity_h(it0,nt,t,x,*args):
    """
    Analytical solution is not derived because of complications.
    Reason: because alpha depends on h.
    Finite difference approximation is used instead.
    """
    D = args[0] # Diffusivity
    K = args[1] # Conductivity
    h = args[2] # Relative boundary conductivity; h = H/K cf. Carslaw and Jaeger, Ch. 1.9 eq. (3)
    l = args[3] # Length of slab
#   Finite difference approx.
    fac = 1.00001
    z = (slab_radiation_unit_production_response_function(it0,nt,t,x,
                                                          *(D,K,h*fac,l)) -
         slab_radiation_unit_production_response_function(it0,nt,t,x,
                                                          *(D,K,h,l)))
    z = z/((fac-1.0)*h)
    return(z)
    
def slab_radiation_unit_production_sensitivity_K(it0,nt,t,x,*args):
    z = np.zeros(nt, dtype=np.float64)
    D = args[0] # Diffusivity
    K = args[1] # Conductivity
    h = args[2] # Relative boundary conductivity; h = H/K cf. Carslaw and Jaeger, Ch. 1.9 eq. (3)
#    l = args[3] # Length of slab
        
    z = (- slab_radiation_unit_production_response_function(it0,nt,t,x,*args)
         + slab_radiation_unit_production_sensitivity_D(it0,nt,t,x,*args)*D
         - slab_radiation_unit_production_sensitivity_h(it0,nt,t,x,*args)*h)
    z = z/K
# #   Finite difference approx.
#     fac = 1.00001
#     z3 = (slab_radiation_unit_production_response_function(it0,nt,t,x,
#                                                 *(D*fac,K*fac,h/fac,l)) -
#          slab_radiation_unit_production_response_function(it0,nt,t,x,
#                                                 *(D,K,h,l)))
#     z3 = z3/((fac-1.0)*K)
#     for it in range(0,100):#nt):
#         if z3[it] != 0.0:
#             print(it,z[it],z3[it],(z[it]-z3[it]),(z[it]-z3[it])/z3[it])
    return(z)  
    
def slab_radiation_unit_production_sensitivity_x(it0,nt,t,x,*args):
    z = np.zeros(nt, dtype=np.float64)
    # if isinstance(args, (tuple,list)):
    D = args[0] # Diffusivity
    K = args[1] # Conductivity
    h = args[2] # Relative boundary conductivity; h = H/K cf. Carslaw and Jaeger, Ch. 1.9 eq. (3)
    l = args[3] # Length of slab
    
#   Use coordinate transform, so z-b.c. wil be at x = 0, not at x = -l, 
#   and zero gradient will be at x = l, not at x = 0.
    xl = x - l

    nmax = 1000
    hl = h*l
    alpha = np.zeros(nmax, dtype=float)
    for n in range(0,nmax):
        alpha[n] = xtan_root(hl, n)/l
        if alpha[n] < 0.0:
            print("\nalpha[%d]=%f < 0.0"%(n,alpha[n]))
    alpha2 = np.power(alpha, 2)
        
    t0 = t[it0]
    h2 = h*h
    Kh = 0.5/(K*h)
    h0 = -2.*h*xl
    EPS = 1.e-12
    for it in range(it0+1,nt):
        at = D*(t[it]-t0)
        z[it] = 1.
        dzsum = 0.0
        dz = 1.
        n = 0
        while np.abs(dz) > EPS and n < nmax:
#            while n < nmax:
            dz = (np.exp(-at*alpha2[n])
                 * (1./(alpha2[n]*((h2+alpha2[n])*l+h)))
                 * (-alpha[n]*np.sin(alpha[n]*xl)/np.cos(alpha[n]*l))
                 )
            dzsum = dzsum + dz
            n += 1
        z[it] = Kh * (h0 - 4.*h2*dzsum)
    # else:
    #     print("*** ERROR *** For slab radiation/production function, parameter values must be given in tupple!")
    return(z)    
    
    
######################################################################################
# Here follows functions to compute response for time series with constant time steps
######################################################################################

# Function to check dimensions on x, t, and bc arrays, and that time steps are constant
#
def cts_pass_check(x,t,bc):

    pass_check = True

    if x.ndim != 1:
        print("*** ERROR *** Dimension of x array cannot exceed 1!")
        pass_check = False

    if t.ndim != 1:
        print("*** ERROR *** Dimension of time array cannot exceed 1!")
        pass_check = False
    else:
        nt = len(t)
        # Check that time steps are constant
        if nt > 1:
            delt = t[1]-t[0]
            fmin = 1.e-10#float_info.epsilon
            for it in range(2,nt):
                if abs(t[it]-t[it-1]-delt) > fmin:
                    print("*** ERROR *** Time steps are not constant!")
                    print(it,(t[it]-t[it-1]),delt,(t[it]-t[it-1])-delt,fmin)
                    pass_check = False

    nbc = bc.shape[0]
    if nbc != nt:
        print("*** ERROR *** Dimension of b.c. array must be the same as dimension of time array!")
        pass_check = False

    return(pass_check)

def plot_response(x,y):
    plt.plot(x,y)
    axs=plt.gca()
    axs.set_xscale("log")
#    axs.set_yscale("log")
    plt.show()
    plt.close()
    
# Function to compute time serie response at nx locations; uses superposition
#
def response(func,bc,nx,x,it0,nt,t,*args,**kwargs):
    z = np.zeros((nt,nx), dtype=np.float64)
    bci = 0.0
    # Initial steady state ?
    for key, value in kwargs.items():
        if str(key).lower() == "bcinit":
            bci = float(value)
            tini = [0.0, 1.e15]
            for ix in range(0,nx):
                zi = func(0,2,tini,x[ix],*args)
                z[:,ix] = bci * zi[1]

#    t1 = time()
#    for ix in range(0,nx):
#        zresp = func(it0-1,nt,t,x[ix],*args)
#        
#        z[:,ix] += (bc[it0-1] - bci) * zresp[:]
#        
#        for it in range(it0,nt):
#            bc0 = bc[it] - bc[it-1]
##            iresp = 1
##            for it1 in range(it,nt): #(it+1,nt):
##                z[it1,ix] += bc0 * zresp[iresp]
##                iresp += 1
#            z[it:nt,ix] += bc0 * zresp[1:nt-it+1] # Use broadcast - MUCH faster!
#    t2 = time()
#    print(ix,t2-t1)

# This approximately halves the computation time compared to the above
    zresp = np.zeros((nt,nx), dtype=np.float64)
    for ix in range(0,nx):
        zresp[:,ix] = func(it0-1,nt,t,x[ix],*args)
    z[:,:] += (bc[it0-1] - bci) * zresp[:,:]
    try:
        if kwargs['aggregate'] == 'sum':
            step_aggregation = False
        else:
            step_aggregation = True
    except:
        step_aggregation = True
    if step_aggregation:
        for it in range(it0,nt):
            bc0 = bc[it] - bc[it-1]
            z[it:nt,:] += bc0 * zresp[1:nt-it+1,:] # Use broadcast - MUCH faster!
    else:
        for it in range(it0,nt):
            z[it:nt,:] += bc[it] * zresp[1:nt-it+1,:] # Use broadcast - MUCH faster!
    return(z)
    
def unit_response_function_dictionary():
    # Returns dictionary of Carslaw and Jaeger unit response functions
    return( {
            "lin_res": lin_res_response,
            "lin_res_sens_tc": lin_res_sens_TC,
        
        
            "solid_step": sinf_solid_unit_step_response,
            "solid_step_sens_d": sinf_solid_unit_step_sensitivity_D,
            "solid_step_sens_x": sinf_solid_unit_step_sensitivity_x,
            
            "solid_lin":  sinf_solid_unity_line_response,
            
            "solid_rad":  sinf_solid_radiation_unit_response_function,
#            "solid_rad_int":  sinf_solid_radiation_unit_response_function_int,
            "solid_rad_sens_d": sinf_solid_radiation_unit_sensitivity_D,
            "solid_rad_sens_h": sinf_solid_radiation_unit_sensitivity_h,
            "solid_rad_sens_x": sinf_solid_radiation_unit_sensitivity_x,

            "solid_prod": sinf_solid_unit_production_response_function,
            "solid_prod_sens_d": sinf_solid_unit_production_sensitivity_D,
            "solid_prod_sens_k": sinf_solid_unit_production_sensitivity_K,
            "solid_prod_sens_x": sinf_solid_unit_production_sensitivity_x,
            
            "solid_rad_prod": sinf_solid_radiation_unit_production_response_function,
            "solid_rad_prod_sens_d": sinf_solid_radiation_unit_production_sensitivity_D,
            "solid_rad_prod_sens_k": sinf_solid_radiation_unit_production_sensitivity_K,
            "solid_rad_prod_sens_h": sinf_solid_radiation_unit_production_sensitivity_h,
            "solid_rad_prod_sens_x": sinf_solid_radiation_unit_production_sensitivity_x,
            
            "slab_step": slab_unit_step_response,
            "slab_step_sens_d": slab_unit_step_sensitivity_D,
            "slab_step_sens_x": slab_unit_step_sensitivity_x,
            
            "slab_prod":  slab_unit_production_response_function,
            "slab_prod_sens_d":  slab_unit_production_sensitivity_D,
            "slab_prod_sens_k":  slab_unit_production_sensitivity_K,
            "slab_prod_sens_x":  slab_unit_production_sensitivity_x,
            
            "slab_rad": slab_radiation_unit_response_function,
            "slab_rad_sens_d": slab_radiation_sensitivity_D,
            "slab_rad_sens_h": slab_radiation_sensitivity_h,
            "slab_rad_sens_l": slab_radiation_sensitivity_l,
            "slab_rad_sens_x": slab_radiation_sensitivity_x,
            
            "slab_rad_prod": slab_radiation_unit_production_response_function,
            "slab_rad_prod_sens_d": slab_radiation_unit_production_sensitivity_D,
            "slab_rad_prod_sens_k": slab_radiation_unit_production_sensitivity_K,
            "slab_rad_prod_sens_h": slab_radiation_unit_production_sensitivity_h,
            "slab_rad_prod_sens_x": slab_radiation_unit_production_sensitivity_x
            })

def response_function_parameters(f):
    par =   {
            "lin_res": tuple(["TC"]),
            "lin_res_sens_tc": tuple(["TC"]),

            "solid_step": tuple(["D"]),
            "solid_step_sens_d": tuple(["D"]),
            "solid_step_sens_x": tuple(["D"]),
            
            "solid_lin":  tuple(["D"]),
            
            "solid_rad":  ("D", "h"),
#            "solid_rad_int":  sinf_solid_radiation_unit_response_function_int,
            "solid_rad_sens_d": ("D", "h"),
            "solid_rad_sens_h": ("D", "h"),
            "solid_rad_sens_x": ("D", "h"),

            "solid_prod": ("D", "K"),
            "solid_prod_sens_d": ("D", "K"),
            "solid_prod_sens_k": ("D", "K"),
            "solid_prod_sens_x": ("D", "K"),
            
            "solid_rad_prod": ("D","K","h"),
            "solid_rad_prod_sens_d": ("D","K","h"),
            "solid_rad_prod_sens_k": ("D","K","h"),
            "solid_rad_prod_sens_h": ("D","K","h"),
            "solid_rad_prod_sens_x": ("D","K","h"),
            
            "slab_step": ("D","l"),
            "slab_step_sens_d": ("D","l"),
            "slab_step_sens_x": ("D","l"),
            
            "slab_prod":  ("D", "K", "l"),
            "slab_prod_sens_d":  ("D", "K", "l"),
            "slab_prod_sens_k":  ("D", "K", "l"),
            "slab_prod_sens_x":  ("D", "K", "l"),
            
            "slab_rad": ("D", "h", "l"),
            "slab_rad_sens_d": ("D", "h", "l"),
            "slab_rad_sens_h": ("D", "h", "l"),
            "slab_rad_sens_l": ("D", "h", "l"),
            "slab_rad_sens_x": ("D", "h", "l"),
            
            "slab_rad_prod": ("D", "K", "h", "l"),
            "slab_rad_prod_sens_d": ("D", "K", "h", "l"),
            "slab_rad_prod_sens_h": ("D", "K", "h", "l"),
            "slab_rad_prod_sens_k": ("D", "K", "h", "l"),
            "slab_rad_prod_sens_x": ("D", "K", "h", "l")
            }
    return(par[f])
    
    
# Function to set unit response function and other variables, then calls response()
#    
def cts_response(bc_type,bc,x,t,*args,**kwargs): # Constant time step response
    
    # Get Dictionary of Carslaw and Jaeger unit response functions
    urfunc = unit_response_function_dictionary()
    bc_type = bc_type.lower()
    nx = len(x)
    nt = len(t)
    it0 = 1
    delt = t[1] - t[0]
    
    if bc_type == "solid_lin":
        # Check for x[ix] == 0.0, which won't work
        for ix in range(0,nx):
            if x[ix] == 0.0:
                print("*** ERROR *** solid linear response cannot be computed for x == 0.0")
                err = True
                return(err)
        # Compute slope for each time step
        bc1 = np.zeros(nt)
        bc1[1] = (bc[1]-bc[0])/delt
        for it in range(2,nt):
            bc1[it] = (bc[it]-bc[it-1])/delt
        bc = bc1
        it0 = 2
        
    try:

        z = response(urfunc[bc_type],bc,nx,x,it0,nt,t,*args,**kwargs)
        
    except KeyError:
        print("\nResponse \"%s\" not implemented!"%bc_type)
        z = np.zeros((nt,nx), dtype=np.float64)

    return(z)
        