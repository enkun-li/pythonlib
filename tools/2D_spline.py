#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: 2D_spline.py
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-05-18 21:00:33
#==================================

import numpy as np

def spline_2D(x, y, xarr, yarr, zmat):
    nx = xarr.shape[0]
    ny = yarr.shape[0]
    
    hx = (xarr[-1] - xarr[0])/(nx -1)
    hy = (yarr[-1] - yarr[0])/(ny -1)
    
    i = int((x-xarr[0])/hx)
    j = int((y-yarr[0])/hy)

    t = (x - xarr[i])/(xarr[i+1] - xarr[i])
    u = (y - yarr[j])/(yarr[j+1] - yarr[j])

    z = (1-t)*(1-u) * zmat[i,j] + t*(1-u) *zmat[i+1,j] \
            + t*u *zmat[i+1,j+1] + (1-t)*u *zmat[i,j+1]
    return z

def bcucof_spline_2D(x, y, xarr, yarr, zmat):
    nx = xarr.shape[0]
    ny = yarr.shape[0]
    
    hx = (xarr[-1] - xarr[0])/(nx -1)
    hy = (yarr[-1] - yarr[0])/(ny -1)
    
    i = int((x-xarr[0])/hx)
    j = int((y-yarr[0])/hy)

    t = (x - xarr[i])/(xarr[i+1] - xarr[i])
    u = (y - yarr[j])/(yarr[j+1] - yarr[j])

    z = (1-t)*(1-u) * zmat[i,j] + t*(1-u) *zmat[i+1,j] \
            + t*u *zmat[i+1,j+1] + (1-t)*u *zmat[i,j+1]
    return z

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ == "__main__":
    xx = np.linspace(0, 2*np.pi, 40)
    yy = np.linspace(0, 2*np.pi, 40)
    zfunc = lambda x,y: x*x+y #np.sin(x) + np.cos(y)
    zz = np.array([[zfunc(x,y) for x in xx] for y in yy])

    print("%12s %12s %12s %12s"%("x", "y", "spline", "real"))
    for i in range(20):
        x = np.random.uniform(0.1,1.99*np.pi)
        y = np.random.uniform(0.1,1.99*np.pi)
        z = spline_2D(x,y,xx,yy,zz)
        print("%12.6f %12.6f %12.6f %12.6f"%(x,y,z,zfunc(x,y)))

