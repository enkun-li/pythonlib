#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: integratefunc.pyx
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-05-17 15:50:28
#==================================

import numpy as np

cdef extern from "math.h":
    double sqrt(double xx)
    double exp(double xx)
    double sinh(double xx)
    double cosh(double xx)
    double sin(double xx)
    double cos(double xx)
    double fabs(double xx)
    double log10(double xx)
    double log(double xx)

cdef:
    int NPOINT = 5
    double PI = 3.14159265359
    double EPS = 1.0e-11

def gauleg(double x1, double x2):
    cdef:
        int m, i, j
        double z1, z, xm, xl, pp, p3, p2, p1
        double[:] x, w
    
    x = np.zeros([NPOINT+1])
    w = np.zeros([NPOINT+1])
    
    m = int( (NPOINT +1)/2 )
    xm = 0.5*(x2+x1)
    xl = 0.5*(x2-x1)
    
    z1 = 1e30
    for i from 1<=i<=m:
        z = cos(PI*(i-0.25)/(NPOINT+0.5))
        pp = 1.0
        while (fabs(z - z1) > EPS):
            p1 = 1.0
            p2 = 0.0
            for j from 1<=j<=NPOINT:
                p3 = p2
                p2 = p1
                p1 = ((2.0*j - 1.0) *z *p2 -(j-1.0) *p3)/j
            pp = NPOINT*(z*p1 -p2)/(z*z -1.0)
            z1 = z
            z = z1 - p1/pp
        x[i] = xm -xl*z
        x[NPOINT+1-i] = xm+xl*z
        w[i] = 2.0*xl/((1.0 -z*z) *pp *pp)
        w[NPOINT+1-i] = w[i]
    return x, w

cpdef double qgausleg(func, double a, double b, args=()):
    cdef:
        double[:] x, w
        double s
        int i

    x, w = gauleg(a, b)
    s = 0
    for i from 1<= i <=NPOINT:
        s += w[i] * func(x[i], *args)
    return s

cpdef double qgausleg2d_fix(func, double a, double b, double c, double d, args=()):
    cdef:
        double[:] xa, xc, wa, wc
        double s
        int i, j

    xa, wa = gauleg(a, b)
    xc, wc = gauleg(c, d)
    s = 0
    for i from 1<= i <=NPOINT:
        for j from 1<= j <= NPOINT:
            s += wa[i] * wc[j] * func(xa[i], xc[j], *args)
    return s

cpdef double qgaus(func, double a, double b, args=()):
    cdef double[:] x, w
    cdef double xm, xr, s, dx
    cdef int i

    if not isinstance(args, tuple):
        args = (args,)

    x = np.array([-0.906179846, -0.538469310, 0.0, 0.538469310, 0.906179846])
    w = np.array([0.236926885, 0.478628670, 0.568888889, 0.478628670, 0.236926885])

    xm = 0.5*(b+a)
    xr = 0.5*(b-a)
    s= 0 
    for i from 0 <= i < 5:
        dx = xr * x[i]
        s += w[i] * func(xm + dx, *args)
    s = s*xr
    return s

#cpdef double dif_qgaus(func, double a, double b, args=()):
#    cdef double[:] xx
#    cdef double s
#    cdef int i, it
#
#    if not isinstance(args, tuple):
#        args = (args,)
#    
#    if(b-a >= 0.5):
#        it = int(b-a)+1
#        xx = np.linspace(a,b,it+1)
#        s = np.sum([qgaus(func, xx[i], xx[i+1], *args) for i in range(it)])
#    else:
#        s = qgaus(func, a, b, *args)
#    return s

cpdef double qgaus2d_fix(func, double a, double b, double c, double d, args=()):
    cdef double[:] x, w
    cdef double xm, xr, dx, ym, yr, dy, s
    cdef int i, j
    
    if not isinstance(args, tuple):
        args = (args,)

    x = np.array([-0.906179846, -0.538469310, 0.0, 0.538469310, 0.906179846])
    w = np.array([0.236926885, 0.478628670, 0.568888889, 0.478628670, 0.236926885])

    xm = 0.5*(b+a)
    xr = 0.5*(b-a)
    ym = 0.5*(d+c)
    yr = 0.5*(d-c)

    s= 0 

    for i from 0 <= i < 5:
        dx = xr * x[i]
        for j from 0 <= j < 5:
            dy = yr * x[j]
            s += w[i] * w[j] * func(xm + dx, ym +dy, *args)
    s = s*xr*yr
    return s

cpdef double qqgaus(func, double a, double b, args=(), double EPS=1e-8, int JMAX=20):
    cdef double s0, s1, err
    cdef int i, j, it
    cdef double[:] xx
    
    if not isinstance(args, tuple):
        args = (args,)

    s0 = qgaus(func, a, b)
    s1 = 0

    for j from 2<=j<JMAX:
        it = 2**j
        xx = np.linspace(a,b,it+1)
        s1 = 0
        for i from 0<=i<it:
            s1 += qgaus(func, xx[i], xx[i+1], *args)
        if(np.abs(s1-s0) < EPS *np.abs(s0)):
            err = np.abs(s1-s0)/np.abs(s0)
            return s1
    print("Too many steps in qqgaus")
    return s1

