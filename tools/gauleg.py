#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: gauleg.py
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-05-17 09:50:34
#==================================

import sys
import numpy as np

EPS = 1e-12

def gauleg(x1, x2, n):
    m = int((n+1)/2)
    xm = 0.5 * (x2 + x1)
    xl = 0.5 * (x2 - x1)

    x = np.zeros(n+1)
    w = np.zeros(n+1)

    z1 = 1e30
    for i in range(1,m+1):
        z = np.cos(np.pi*(i-0.25)/(n+0.5))
        while (np.abs(z - z1) > EPS ):
            p1 = 1.0
            p2 = 0.0
            for j in range(1,n+1):
                p3 = p2
                p2 = p1
                p1 = ((2.0*j - 1.0)*z*p2 -(j-1.0)*p3)/j
            pp = n*(z*p1 -p2)/(z*z-1.0)
            z1 = z
            z = z1 - p1/pp
        x[i] = xm -xl*z
        x[n+1-i] = xm +xl*z
        w[i] = 2.0*xl/((1.0-z*z)*pp*pp)
        w[n+1-i] = w[i]
    return x, w

def qgauleg(func, a, b, n=5):
    x, w = gauleg(a,b,n)
    ss = 0
    for i in range(1,n+1):
        ss += w[i] * func(x[i])
    return ss

def qqgauleg(func, a, b, n=5):
    s0 = qgauleg(func, a, b, n)
    s1 = qgauleg(func, a, b, 2*n)
    print("="*60)
    print("n= %3d value is %12.6f"%(n, s0))
    print("n= %3d value is %12.6f"%(2*n, s1))
    print("The difference between n= %d and 2*n = %d "%(n, 2*n), end=' ')
    print("is %20.16f"%(s1-s0))
    return 

def dqgauleg(func, a, b, n=5):
    s0 = qgauleg(func, a, b, n)
    JMAX = 20
    EPS = 1e-7
    
    rs = actual_func(a,b)

    for j in range(2,JMAX):
        #h = (b-a)/2**j
        it = 2**j
        xx = np.linspace(a,b,it+1)
        s1 = 0
        for i in range(it):
            #x1 = a + i*h
            #x2 = a + (i+1)*h
            s1 += qgauleg(func, xx[i], xx[i+1], n)
        
        print("="*60)
        print("j= %3d value is %12.6f"%(j-1, s0))
        print("j= %3d value is %12.6f"%(j, s1))
        print("The difference is %20.16e"%((s1-s0)/s0))
        print("The err is        %20.16e"%((s1-rs)/rs))

        if (np.abs(s1-s0) < EPS*np.abs(s0)):
            return s1
        s0 = s1
    print("Can not get the precision!")
    return s1

def arr_gauleg(func, a, b, n=5):
    JMAX = int((b-a))+1
    print("="*60)
    print("JMAX = %d"%JMAX)
    #h = (b-a)/JMAX
    xx = np.linspace(a,b,JMAX+2)
    
    ss = 0
    for j in range(JMAX):
        ss += qgauleg(func, xx[j], xx[j+1], n)
    
    return ss


def print_x_w(x,w,n):
    print("%5s %12s %12s"%('#', 'x[i]', 'w[i]'))
    for i in range(n+1):
        print('%5d %12.9f %12.9f'%(i, x[i], w[i]))
    return

def func(x):
    return x*np.exp(-x)

def actual_func(x1,x2):
    return (1+x1)*np.exp(-x1) -(1+x2)*np.exp(-x2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__=="__main__":
    arg = sys.argv[1:]
    x1 = float(arg[0])
    x2 = float(arg[1])
    n = int(arg[2])
    
    x, w = gauleg(x1,x2,n)
    print_x_w(x,w,n)

    ss = qgauleg(func, x1, x2, n)
    rs = actual_func(x1, x2)
    print("="*60)
    print("Actual value is %12.6f"%rs)
    print("Integrate value %12.6f"%ss)
    qqgauleg(func, x1, x2, n)
    
    dqgauleg(func, x1, x2, n)

    srs = arr_gauleg(func, x1, x2, n)
    print("Integrate value %12.6f"%srs)
    print("err with actual %12.6e"%((srs-rs)/rs))
