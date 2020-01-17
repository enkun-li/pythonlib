#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: mcmc.pyx
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-01-13 19:29:15
#==================================

cimport cython
import numpy as np
import scipy.optimize as opt
from scipy import integrate

cdef extern from "math.h":
    double sqrt(double xx)
    double exp(double xx)
    double sinh(double xx)
    double sin(double xx)
    double cos(double xx)
    double fabs(double xx)
    double log10(double xx)
    double log(double xx)
#
# mcmc
#
def mcmc(aimfunc, double[:,:] initheta, int nsample=500):
    '''
    '''
    cdef double[:,:] result
    cdef double[:] pre_theta, pro_theta
    cdef double[:,:] cov
    cdef double pre_aim, pro_aim, alpha, uni
    cdef int i, npar
    cdef double Ti, Tmax, Tmin, ratio

    npar = np.shape(initheta)[0]
    pre_theta = np.random.uniform(initheta[:,1], initheta[:,2])
    cov = np.eye(npar)*initheta[:,3]
    pre_aim = aimfunc(pre_theta,initheta[:,1],initheta[:,2])

    result = np.zeros((npar+1, nsample))
    
    Tmax = 100
    Tmin = 1e-5
    ratio = 0.68

    Ti = Tmax
    while(Ti > Tmin):
        for i from 0<= i <nsample:
            pro_theta = np.random.multivariate_normal(pre_theta, cov)
            pro_aim = aimfunc(pro_theta, initheta[:,1], initheta[:,2])
            
            if(pro_aim > pre_aim):
                pre_theta = pro_theta
                pre_aim = pro_aim
            else:
                #alpha = np.min([1., np.exp(pro_aim -pre_aim)] )
                alpha = 1.0/(1.0+exp(-(pro_aim -pre_aim)/Ti) )
                uni = np.random.uniform()
                if( uni < alpha):
                    pre_theta = pro_theta
                    pre_aim = pro_aim
                    
            result[0,i] = pre_aim
            result[1:,i] = pre_theta[:]
        Ti = Ti*ratio
    return result

def aimfunc(double[:] theta, double[:] sigdn, double[:] sigup):
    cdef double chisq, sq, ns
    cdef int n

    n = np.shape(theta)[0]
    chisq = 0.0
    for i from 0 <= i <n:
        if(sigdn[i] < theta[i] < sigup[i]):
            chisq += theta[i]**2 - np.cos(10*theta[i])
        else:
            chisq = 1e30
    return chisq
