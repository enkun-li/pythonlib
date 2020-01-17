#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: qgaus.py
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-05-17 15:16:33
#==================================

import numpy as np
import sys
from scipy.integrate import quad, dblquad
import time
from integratefunc import qqgaus as int_qqgaus
from integratefunc import qgaus as int_qgaus
from integratefunc import qgaus2d_fix

def qgaus(func, a, b):
    x = [-0.906179846, -0.538469310, 0.0, 0.538469310, 0.906179846]
    w = [0.236926885, 0.478628670, 0.568888889, 0.478628670, 0.236926885]

    xm = 0.5*(b+a)
    xr = 0.5*(b-a)
    s= 0 
    for i in range(5):
        dx = xr * x[i]
        s += w[i] * func(xm + dx)
    return s*xr

def qqgaus(func, a, b):
    s0 = qgaus(func, a, b)
    JMAX = 20
    EPS = 1e-8

    for j in range(2,JMAX):
        it = 2**j
        xx = np.linspace(a,b,it+1)
        s1 = 0
        for i in range(it):
            s1 += qgaus(func, xx[i], xx[i+1])
        if(np.abs(s1-s0) < EPS *np.abs(s0)):
            err = np.abs(s1-s0)/np.abs(s0)
            return s1, err
    print("Too many steps in qqgaus")
    return s1, 0.0


def func(x):
    return np.sqrt(1+x*x)

def funcxy(x,y):
    ff = np.log(x+2*y)
    return ff

def funcxy_args(x,y, alpha, beta):
    ff = np.log(alpha*x+beta*y)
    return ff

def real_int_func(x):
    f = 0.5*(x*np.sqrt(1+x*x) +np.arcsinh(x))
    return f

#***********************#
if __name__ == "__main__":
    arg = sys.argv[1:] 
    a = float(arg[0])
    b = float(arg[1])

    rss = real_int_func(b) -real_int_func(a)
    print("Real result is %20.15f\n"%(rss))

    ts = time.time()
    ss = qgaus(func, a, b)
    te = time.time()
    print("qgaus result is      %20.15f"%ss)
    print("Precision is         %20.15e"%((ss-rss)/rss))
    print("Cost is %10.4e\n"%(te - ts))

    ts = time.time()
    ss = int_qgaus(func, a, b)
    te = time.time()
    print("int_qgaus result is  %20.15f"%ss)
    print("Precision is         %20.15e"%((ss-rss)/rss))
    print("Cost is %10.4e\n"%(te - ts))

    ts = time.time()
    ss, err = quad(func, a, b)
    te = time.time()
    print("quad result is       %20.15f"%ss)
    print("Precision is         %20.15e"%((ss-rss)/rss))
    print("Cost is %10.4e\n"%(te - ts))
    
    ts = time.time()
    ss, err = qqgaus(func, a, b)
    te = time.time()
    print("qqgaus result        %20.15f"%ss)
    print("Precision is         %20.15e"%((ss-rss)/rss))
    print("Cost is %10.4e\n"%(te - ts))
    
    ts = time.time()
    ss = int_qqgaus(func, a, b)
    te = time.time()
    print("int_qqgaus result    %20.15f"%ss)
    print("Precision is         %20.15e"%((ss-rss)/rss))
    print("Cost is %10.4e\n"%(te - ts))

    
    #===========================================
    ts = time.time()
    ss = qgaus2d_fix(funcxy, 1.4, 2.0, 1.0, 1.5)
    te = time.time()
    print("qgaus2d_fix result   %20.15f"%ss)
    print("Cost is %10.4e\n"%(te - ts))
    
    ff = lambda y,x: funcxy(x,y)

    ts = time.time()
    ss, err = dblquad(ff, 1.4, 2.0, 1.0, 1.5)
    te = time.time()
    print("dblquad result       %20.15f"%ss)
    print("Cost is %10.4e\n"%(te - ts))
    
    #****************************
    ts = time.time()
    ss = qgaus2d_fix(funcxy_args, 1.4, 2.0, 1.0, 1.5, args=(1,2))
    te = time.time()
    print("qgaus2d_fix result   %20.15f"%ss)
    print("Cost is %10.4e\n"%(te - ts))

