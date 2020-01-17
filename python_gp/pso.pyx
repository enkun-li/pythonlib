#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: pso.pyx
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-01-13 19:28:20
#==================================

import numpy as np

#********************************************************
# particle swarm optimization algorithm:
#
cpdef double[:,:] pso(lnlike, double[:,:] dats, double[:] sigma, \
        int nsample=500, int npar=2):
    '''
    cpdef double[:,:] pso(lnlike, double[:,:] dats, double[:] sigma, int nsample=500, int npar=2):
    ====================

    '''
    cdef double[:,:] Xt, Vt, result
    cdef double[:,:] Pbest
    cdef double[:] P_fit, Gbest
    cdef double G_fit, tmp
    cdef int npart, i, j, k
    cdef double[:] Xtt, Vtt, xmin, xmax, vmax
    cdef double r1, r2, c1, c2, w, wmax, wmin, t

    npart = 50
    result = np.zeros((2, npar))
    Xt = np.zeros((npart, npar))
    Vt = np.zeros((npart, npar))
    Xtt = np.zeros(npar)
    Vtt = np.zeros(npar)
    xmax = np.zeros(npar)
    xmin = np.zeros(npar)
    vmax = np.zeros(npar)

    Pbest = np.zeros((npart, npar))
    Gbest = np.zeros(npar)
    P_fit = np.zeros(npart)
    G_fit = 1e30
    
    for i from 0<= i < npar:
        xmax[i] = sigma[i]
        xmin[i] = 0.01
        vmax[i] = 0.5*(xmax[i] -xmin[i])

    wmax = 0.9
    wmin = 0.1
    t = 0.0
    c1 = 1.193
    c2 = 1.193
    w = wmax

    for i from 0<= i < npart:
        for j from 0<= j < npar:
            Xt[i,j] = np.random.uniform(xmin[j], xmax[j])
            Vt[i,j] = np.random.uniform(-1.0,1.0)*vmax[j]
            Pbest[i,j] = Xt[i,j]
        P_fit[i] = -lnlike(Pbest[i,:], sigma, dats)
        if(P_fit[i] < G_fit):
            G_fit = P_fit[i]
            Gbest[:] = Pbest[i,:]

    while w >= wmin:
        for i from 0<= i < nsample:
            for j from 0<= j< npart:
                r1 = np.random.uniform()
                r2 = np.random.uniform()
                for k from 0<= k < npar:
                    Xtt[k] = Xt[j,k] + Vt[j,k]
                    Vtt[k] = w*Vt[j,k] +c1*r1*(Pbest[j,k]-Xt[j,k]) \
                            +c2*r2*(Gbest[k]-Xt[j,k])
                    if(Xtt[k] < xmin[k]):
                        Xt[j,k] = xmin[k]
                        Vt[j,k] = -Vtt[k]
                    elif(Xtt[k] > xmax[k]):
                        Xt[j,k] = xmax[k]
                        Vt[j,k] = -Vtt[k]
                    elif(Vtt[k] > vmax[k]):
                        Vt[j,k] = vmax[k]
                    elif(Vtt[k] < -vmax[k]):
                        Vt[j,k] = -vmax[k]
                    else:
                        Xt[j,k] = Xtt[k]
                        Vt[j,k] = Vtt[k]
                tmp = -lnlike(Xt[j,:], sigma, dats)
                if(tmp < P_fit[j]):
                    P_fit[j] = tmp
                    Pbest[j,:] = Xt[j,:]
                    if(P_fit[j] < G_fit):
                        G_fit = P_fit[j]
                        Gbest[:] = Pbest[j,:]
        t += 1
        w = wmax - (t-1)/10.0*(wmax -wmin)
        print('='*50)
        print('The weight parameter is: %8.2f'%w)
        print('The best fit is: %9.4f'%G_fit)
        print('The best position is: ')
        for i from 0<=i<npar:
            print('%s-th:  %9.4f'%(i, Gbest[i]))
        print('')

    result[0,0] = G_fit
    result[1,:] = Gbest[:]
    return result
