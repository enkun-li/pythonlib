#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: special_functions.pyx
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-03-27 08:52:08
#==================================

import numpy as np

def GaussianDistribution(x, mu, sig ):
    '''
    multi_GaussianDistribution(x,mu,cov)
    =====================================
    x: single or array like values
    mu: single for 1-D
    sig: single for 1-D
    '''    
    p = lambda y: 1/np.sqrt(2*np.pi*sig**2)*np.exp( -0.5*(y-mu)**2/sig**2)
    pdf = np.array([p(y) for y in x])
    return pdf

def multi_GaussianDistribution(x, double[:] mu, double[:,:] cov ):
    '''
    multi_GaussianDistribution(x,mu,cov)
    =====================================
    x:  array like values
    mu: array like for N-D
    cov: matrix for N-D
    '''
    cdef double[:,:] L, B, Eye, M
    cdef double[:] vec
    cdef double logdetM, detM, pdf
    n = np.shape(mu)[0]
    Eye = np.ones(n)
    L = np.linalg.cholesky(cov)
    B = np.linalg.solve(L,Eye)
    M = np.linalg.solve(np.transpose(L), B)
    logdetM = np.sum(np.log(np.diagonal(L))) # *2
    detM = np.exp(logdetM)
    vec = np.array([xx-mm for xx, mm in zip(x,mu)])
    chisq = np.dot(np.matmul(vec, M), vec)
    pdf = np.exp(-0.5*chisq)/(2*np.pi)**(n/2.0)/detM
    return pdf

cpdef double cubic_spline(double x, double[:] xx, double[:] yy,double[:] yy2,
        int order=0):
    '''
    cubic_spline(x, xx, yy, yy2)
    =====================================
    From array {x,f(x),f''(x)} to generate f(z)
    param x: the aim position, steep must be equal!
    param xx: x axis of data
    param yy: f(x) of data
    param yy2: second derivation of f(x) => f''(x)
    param order: 0 for f(x)
                 1 for f'(x)
                 2 for f''(x)
    '''
    cdef int n, j
    cdef double h, A, B
    n = xx.shape[1]
    h = (xx[-1]-xx[0])/(n-1)
    j = int((x-xx[0])/h)
    B = (x-xx[j])/h
    A = (xx[j+1]-x)/h
    if(order == 0):
        fx = A*yy[j] +B*yy[j+1] +h**2/6*( (A**3-A)*yy2[j] +(B**3-B)*yy2[j+1] )
    elif(order==1):
        fx = (yy[j+1]-yy[j])/h -(3*A**2-1)/6*h*yy2[j] +(3*B**2-1)/6*h*yy2[j+1]
    elif(order==2):
        fx = A*yy2[j] +B*yy2[j+1]
    else:
        print('Can not calculate %d order.'%order)
    #print(n, h, j)
    return fx

cpdef double[:] spline(double[:] xx, double[:] yy, double yp0=0, double ypn=0):
    '''
    spline(xx,yy,yp0,ypn)
    =============================
    generate the second dervation of array.
    '''
    cdef int n, k, i
    cdef double[:] y2, u
    cdef double sig, p, qn, un

    n = xx.shape[0]
    y2 = np.zeros(n)
    u = np.zeros(n+1)
    if(yp0>1e30):
        y2[0]=0
        u[0] = 0
    else:
        y2[0] = -0.5
        u[0] = (3/(xx[1] -xx[0]))* ((yy[1]-yy[0])/(xx[1]-xx[0]) -yp0)
    for i in range(1,n-1):
        sig = (xx[i]-xx[i-1])/(xx[i+1]-xx[i-1])
        p = sig*y2[i-1] +2
        y2[i] = (sig-1)/p
        u[i] = (6.*((yy[i+1]-yy[i])/(xx[i+1]-xx[i])-(yy[i]-yy[i-1])/
            (xx[i]-xx[i-1]))/(xx[i+1]-xx[i-1])-sig*u[i-1])/p
    if(ypn > 1e30):
        qn =0
        un = 0
    else:
        qn = 0.5
        un = (3/(xx[n-1] -xx[n-2]))*(ypn-(yy[n-1]-yy[n-2])/(xx[n-1]-xx[n-2]))
        y2[n-1] = (un-qn*u[n-2])/(qn*y2[n-2]+1)
    for k in np.arange(n-2,0,-1):
        y2[k] = y2[k]*y2[k+1] +u[k]        
    return y2
