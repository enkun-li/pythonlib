#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: gapp_kernel.pyx
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2018-12-25 21:29:48
#==================================

cimport cython

import numpy as np
import scipy.optimize as opt
from scipy import integrate
from gapp import dgp

#@cython.boundscheck(False)
#@cython.wraparound(False)

cdef extern from "math.h":
    double sqrt(double xx)
    double exp(double xx)
    double sinh(double xx)
    double sin(double xx)
    double cos(double xx)
    double fabs(double xx)
    double log10(double xx)
    double log(double xx)

def load_dat(filename):
    '''
    read in data:
        load_dat(filename)
    filename: the name of data file
    '''
    cdef double[:,:] dats
    dats = np.loadtxt(filename, unpack=True)
    return dats


cpdef double kernel(double x, double y, double sigmaf=1, \
        double length=2):
    '''
    Gaussian kernel: sig**2 exp(-(x-y)**2/2/l**2)
    '''
    return sigmaf**2*exp(-(x-y)**2/2/length**2)

cpdef double dkerneldx(double x, double y, double sigmaf=1, \
        double length=2):
    '''
    dk(x,y)/dx
    '''
    return (y-x)/length**2*kernel(x,y,sigmaf,length)

cpdef double dkerneldy(double x, double y, double sigmaf=1,\
        double length=2):
    '''
    dk(x,y)/dy
    '''
    return (x-y)/length**2*kernel(x,y,sigmaf,length)

cpdef double d2kerneldxdy(double x, double y, double sigmaf=1, \
        double length=2):
    '''
    d^2k(x,y)/dxdy
    '''
    return (1-(x-y)**2/length**2)/length**2*kernel(x,y,sigmaf,length)

cpdef double[:,:] Mcov(double[:,:] dats, double sigmaf=1, double length=2):
    '''
    M: K(X,X)+C
        Mcov(dats,sigf,lenf)
    dats: the datas
    '''
    cdef double[:,:] cov
    cdef double diag
    cdef int i, j, n
    n = np.shape(dats)[1]
    cov = np.zeros((n,n))
    for i from 0<= i < n:
        for j from 0<= j < n:
            cov[i,j] = kernel(dats[0,i],dats[0,j],sigmaf,length)
        diag = cov[i,i] + dats[2,i]**2.0
        cov[i,i] = diag
    return cov

cpdef double[:,:] InvMat(double[:,:] dats, double sigmaf=1, double length=2):
    '''
    Inv(M): M = K(X,X) +C
    return M & M^{-1}
    '''
    cdef double[:,:] cov, invcov

    #num = np.shape(datx)[0]
    cov = Mcov(dats,sigmaf,length)
    invcov = np.linalg.inv(cov)
    return invcov

cpdef double meanf(double x, double[:,:] dats, double[:,:] invcov,\
        double sigmaf=1, length=2):
    '''
    mean: x = x^*, dats=dats, invcov = InvMat
        meanf(x,dats,invcov,sigmaf,length)
    return mean
    '''
    cdef int num, i
    cdef double[:] vec

    num = np.shape(dats)[1]
    vec = np.zeros(num)
    for i from 0<= i <num:
        vec[i] = kernel(x, dats[0,i], sigmaf, length)    
    return np.dot( np.matmul(vec, invcov), dats[1,:]) 

cpdef double dmeanf(double x, double[:,:] dats, double[:,:] invcov, \
        double sigmaf=1, double length=2):
    cdef int num
    cdef double[:] vec

    num = np.shape(dats)[1]
    vec = np.zeros(num)
    for i in range(num):        
        vec[i] = dkerneldx(x, dats[0,i], sigmaf, length)

    return np.dot( np.matmul(vec, invcov), dats[1,:]) 

cpdef double covariance(double x, double y, double[:,:] dats, \
        double[:,:] invcov, double sigmaf=1, double length=2):
    cdef int num
    cdef double[:] vec1, vec2
    num = np.shape(dats)[1]
    vec1 = np.zeros(num)
    vec2 = np.zeros(num)
    for i from 0<= i <num:
        vec1[i] = kernel(x,dats[0,i],sigmaf, length)
        vec2[i] = kernel(dats[0,i], y, sigmaf, length)
    return kernel(x,y,sigmaf,length) - np.dot(np.matmul(vec1,invcov),vec2)

cpdef double sigma(double x, double[:,:] dats, double[:,:] invcov, \
        double sigmaf=1, double length=2):
    return np.sqrt(covariance(x,x,dats, invcov, sigmaf,length))

cpdef double covariance_dfdf(double x, double y, double[:,:] dats, \
        double[:,:] invcov, double sigmaf=1, double length=2):
    '''
    covariance_dfdf: cov(df,df)
        covariance_dfdf(x,dats,invcov,sigmaf,length)
        the covariance between (df(x), df(y))
    x: the position
    dats: the data
    invcov: Inv(M)
    sigmaf:
    length:
    '''
    cdef int num, i
    cdef double[:] vec1, vec2
    cdef double sig2

    num = np.shape(dats)[1]
    vec1 = np.zeros(num)
    vec2 = np.zeros(num)
    for i from 0<= i <num:
        vec1[i] = dkerneldx(x, dats[0,i], sigmaf, length)
        vec2[i] = dkerneldy(dats[0,i], y, sigmaf, length)
    sig2 = d2kerneldxdy(x,y,sigmaf,length) -np.dot(np.matmul(vec1, invcov), vec2)
    return sig2

cpdef double covariance_fdf(double x, double y, double[:,:] dats, \
        double[:,:] invcov, double sigmaf=1, double length=2):
    '''
    covariancefdf: cov(f,df)
        covariancefdf(x,dats,invcov,sigmaf,length)
        the covariance between (f(x), df(y))
    '''
    cdef int num, i
    cdef double[:] vec1, vec2
    cdef double sig2
    num = np.shape(dats)[1]
    vec1 = np.zeros(num)
    vec2 = np.zeros(num)
    for i from 0<= i <num:
        vec1[i] = kernel(x, dats[0,i], sigmaf, length)
        vec2[i] = dkerneldx(x, dats[0,i], sigmaf, length)
    sig2 = dkerneldy(x,x, sigmaf, length) - np.dot( np.matmul(vec1, invcov), vec2)
    return sig2

cpdef double[:] arr_meanf(double[:] arrx, double[:,:] dats, double[:,:] invcov, \
        double sigmaf=1, double length=2, int default=0):
    '''
    arr_meanf: array of meanf
        arr_meanf(arrx, dats,invcov,sigmaf,length,default)
    arrx: array of x
    default: 0 => mean(f)
             1 => mean(df)
    '''
    cdef int num
    cdef double[:] arrf
    num = np.shape(arrx)[0]
    arrf = np.zeros(num)
    for i from 0<= i <num:
        if(default==0):
            arrf[i] = meanf(arrx[i], dats, invcov, sigmaf, length)
        elif(default==1):
            arrf[i] = dmeanf(arrx[i], dats, invcov, sigmaf, length)
    return arrf

cpdef double[:] arr_sigma(double[:] arrx, double[:,:] dats, double[:,:] invcov, \
        double sigmaf=1, double length=2, int default=0):
    '''
    arr_sigma: array of sigmaf
        arr_sigma(arrx,dats,invcov,sigmaf,length,default)
    arrx: array of x
    default: 0 => sigma(f)
             1 => sigma(df)
    '''
    cdef int num
    cdef double[:] arrf
    num = np.shape(arrx)[0]
    arrf = np.zeros(num)
    for i from 0<= i <num:
        if(default==0):
            arrf[i] = sigma(arrx[i], dats, invcov, sigmaf, length)
        elif(default==1):
            arrf[i] = np.sqrt(covariance_dfdf(arrx[i], arrx[i], dats,invcov,sigmaf, length) )
    return arrf

#**************************************
# likelihood
#
cpdef double lnlikelihood(double[:] theta, double[:,:] dats):
    '''
    lnlikelihood(double[:] theta, double[:,:] dats):
    cdef double[:,:] cov
    cdef double[:,:] invcov
    cdef double sigmaf, length, detcov, lnlike
    cdef int num
    sigmaf, length = theta
    num = np.shape(dats)[0]
    cov = Mcov(dats, sigmaf, length)
    invcov = InvMat(dats,sigmaf,length)
    detcov = np.linalg.det(cov)
    #theory = arr_meanf(datx, datx, daty, date, sigmaf, length, default)
    #lnlike = -0.5*np.dot(np.matmul(daty-theory, invcov), daty-theory) \
    #    - 0.5*np.log(np.abs(detcov)) -num/2*np.log(2*np.pi)
    lnlike = -0.5*np.dot(np.matmul(dats[1,:], invcov), dats[1,:]) \
        - 0.5*log(np.abs(detcov)) -num/2*log(2*np.pi)
    return lnlike
    '''
    cdef double[:,:] cov
    cdef double[:,:] invcov
    cdef double sigmaf, length, detcov, lnlike
    cdef int num
    sigmaf, length = theta
    num = np.shape(dats)[0]
    cov = Mcov(dats, sigmaf, length)
    invcov = InvMat(dats,sigmaf,length)
    detcov = np.linalg.det(cov)
    #theory = arr_meanf(datx, datx, daty, date, sigmaf, length, default)
    #lnlike = -0.5*np.dot(np.matmul(daty-theory, invcov), daty-theory) \
    #    - 0.5*np.log(np.abs(detcov)) -num/2*np.log(2*np.pi)
    lnlike = -0.5*np.dot(np.matmul(dats[1,:], invcov), dats[1,:]) \
        - 0.5*log(np.abs(detcov)) -num/2*log(2*np.pi)
    return lnlike

cpdef double lnprior(double[:] theta, double[:] sigma):
    '''
    lnprior(double[:] theta, double[:] sigma):
    cdef double sigf, lenf
    sigf, lenf = theta
    if(0<sigf<sigma[0] and 0<lenf<sigma[1]):
        return 0.0
    return -1.0e30
    -np.inf
    '''
    cdef double sigf, lenf
    sigf, lenf = theta
    if(0<sigf<sigma[0] and 0<lenf<sigma[1]):
        return 0.0
    return -1.0e30

#***************************************
# lnprobability
#
cpdef double lnprobability(double[:] theta, double[:] sigma, double[:,:] dats):
    '''
    lnprobability(double[:] theta, double[:] sigma, double[:,:] dats):
    cdef double lp
    lp = lnprior(theta,sigma)
    if(not np.isfinite(lp)):
        return -np.inf
    return lnlikelihood(theta, dats)+lp
    '''
    cdef double lp
    lp = lnprior(theta,sigma)
    if(not np.isfinite(lp)):
        return -np.inf
    return lnlikelihood(theta, dats)+lp

#****************************************
# gapp to give sigmaf and length
#
cpdef double[:] output_gp_pars(double[:,:] dats):
    '''
    cpdef double[:] output_gp_pars(double[:,:] dats):
    =================================================
    cdef double zi, zx
    cdef int ns
    cdef double[:] initheta
    cdef double[:] theta
    cdef double[:,:] recH

    theta = np.zeros(2,dtype=np.float64)

    ns = 2
    zi = dats[0,0]
    zx = dats[0,-1]
    recH = np.zeros((ns,3),dtype=np.float64)

    initheta = np.array([100,2], dtype=np.float64)
    gpr = dgp.DGaussianProcess(dats[0],dats[1],dats[2],cXstar=(zi,zx,ns),\
            theta=initheta)
    recH, theta = gpr.gp()
    return theta
    '''
    cdef double zi, zx
    cdef int ns
    cdef double[:] initheta
    cdef double[:] theta
    cdef double[:,:] recH

    theta = np.zeros(2,dtype=np.float64)

    ns = 2
    zi = dats[0,0]
    zx = dats[0,-1]
    recH = np.zeros((ns,3),dtype=np.float64)

    initheta = np.array([100,2], dtype=np.float64)
    gpr = dgp.DGaussianProcess(dats[0],dats[1],dats[2],cXstar=(zi,zx,ns),\
            theta=initheta)
    recH, theta = gpr.gp()
    return theta

#**************************
# rec Hz
cpdef double[:] rec_Hz(double[:] z, double[:,:] dats, double[:,:] invcov, \
        double sigf, double lenf):
    '''
    cpdef double[:] Hz(double[:] z, double[:,:] dats, double[:,:] invcov,
        double sigf, double lenf):
    =================================================================
    cdef double[:] arr_Hz
    cdef int n, i
    
    n = np.shape(z)[0]
    arr_Hz = np.zeros(n,dtype=np.float64)

    for i from 0<= i < n:
        arr_Hz[i] = meanf(z, dats, invcov, sigf, lenf)

    return arr_Hz
    '''
    cdef double[:] arr_Hz
    cdef int n, i
    
    n = np.shape(z)[0]
    arr_Hz = np.zeros(n,dtype=np.float64)

    for i from 0<= i < n:
        arr_Hz[i] = meanf(z[i], dats, invcov, sigf, lenf)

    return arr_Hz

#*****************************
# cov(Hx,Hy)
#
cpdef double[:,:] cov_Hz(double[:] arr_z, double[:,:] dats, \
        double[:,:] invcov, double sigf, double lenf):
    '''
    cpdef double[:,:] cov_Hz(double[:] arr_z, double[:,:] dats, 
        double[:,:] invcov, double sigf, double lenf)
    ==========================================================
    cdef int i,j,n
    cdef double[:,:] covhz

    n = np.shape(arr_z)[0]
    covhz = np.zeros((n,n), dtype=np.float64)
    for i from 0 <= i < n:
        for j from 0 <= j < n:
            covhz[i,j] = covariance(arr_z[i],arr_z[j],dats,invcov,sigf,lenf)
            
    return covhz
    '''
    cdef int i,j,n
    cdef double[:,:] covhz

    n = np.shape(arr_z)[0]
    covhz = np.zeros((n,n), dtype=np.float64)
    for i from 0 <= i < n:
        for j from 0 <= j < n:
            covhz[i,j] = covariance(arr_z[i],arr_z[j],dats,invcov,sigf,lenf)
            
    return covhz

#***************************
# rec dc
#
cpdef double invHz(double z, double[:,:] dats, double[:,:] invcov, \
        double sigf, double lenf):
    '''
    cpdef double invHz(double z, double[:,:] dats, double[:,:] invcov, 
        double sigf, double lenf)
    ===========================================================
    return 1.0/meanf(z,dats,invcov,sigf,lenf)
    '''
    return 1.0/meanf(z,dats,invcov,sigf,lenf)

cpdef double[:] rec_dc(double[:] z, double[:,:] dats, double[:,:] invcov, \
        double sigf, double lenf):
    '''
    cpdef double[:] rec_dc(double[:] z, double[:,:] dats, double[:,:] invcov,
        double sigf, double lenf)
    =============================================
    cdef int n, i
    cdef double[:] recDc
    cdef double fz, err

    n = np.shape(z)[0]
    recDc = np.zeros(n,dtype=np.float64)

    for i from 0 < i< n:
        recDc[i], err = integrate.quad(invHz, 0.0, z[i], \
                args=(dats,invcov,sigf,lenf))

    return recDc
    '''
    cdef int n, i
    cdef double[:] recDc
    cdef double fz, err

    n = np.shape(z)[0]
    recDc = np.zeros(n,dtype=np.float64)

    for i from 0 < i< n:
        recDc[i], err = integrate.quad(invHz, 0.0, z[i], \
                args=(dats,invcov,sigf,lenf))

    return recDc


#****************************************
# \int 1/Hx^2 1/Hy^2 cov(Hx,Hy) dx dy
#
cpdef double int_HxHy(double x, double y, double[:,:] dats, \
        double[:,:] invcov, double sigf, double lenf):
    '''
    cpdef double int_HxHy(double x, double y, double[:,:] dats, 
        double[:,:] invcov, double sigf, double lenf)
    ==========================================================
    cdef double covhxhy, hx, hy

    covhxhy = covariance(x,y,dats,invcov,sigf,lenf)
    hx = meanf(x,dats,invcov,sigf,lenf)
    hy = meanf(y,dats,invcov,sigf,lenf)
    return 1.0/hx**2/hy**2*covhxhy
    '''
    cdef double covhxhy, hx, hy
    cdef double inthxhy

    covhxhy = covariance(x,y,dats,invcov,sigf,lenf)
    hx = meanf(x,dats,invcov,sigf,lenf)
    hy = meanf(y,dats,invcov,sigf,lenf)
    inthxhy = 1.0/hx**2/hy**2*covhxhy
    return inthxhy

cpdef double cov_dcij(double zi, double zj, double[:,:] dats, \
        double[:,:] invcov, double sigf, double lenf):
    '''
    '''
    cdef double covij, err

    covij, err = integrate.dblquad(int_HxHy, 0,zi, 0, zj, \
            args=(dats,invcov,sigf,lenf))
    return covij

cpdef double[:,:] cov_dc(double[:] arr_z, double[:,:] dats, \
        double[:,:] invcov, double sigf, double lenf):
    '''
    cpdef double[:,:] cov_dc(double[:] arr_z, double[:,:] dats, 
        double[:,:] invcov, double sigf, double lenf)
    ==========================================================
    cdef int i,j, n
    cdef double[:,:] covdc

    n = np.shape(arr_z)
    covdc = np.zeros((n,n),dtype=np.float64)

    for i from 0<= i < n:
        for j from 0<= j <n:
            covdc[i,j] = cov_dcij(arr_z[i],arr_z[j],dats,invcov,sigf,lenf)

    return covdc
    '''
    cdef int i,j, n
    cdef double[:,:] covdc

    n = np.shape(arr_z)[0]
    covdc = np.zeros((n,n),dtype=np.float64)

    for i from 0<= i < n:
        for j from 0<= j <n:
            covdc[i,j] = cov_dcij(arr_z[i],arr_z[j],dats,invcov,sigf,lenf)

    return covdc


#********************************************************
# mcmc
#
cpdef double[:,:] mcmc(lnlike, double[:,:] dats, double[:] initheta, \
        double[:] sigma, int nsample=500, int npar=2):
    '''
    mcmc(lnlike, double[:,:] dats, double[:] initheta,
        double[:] sigma, int nsample=500, int npar=2):
    cdef double[:,:] result
    cdef double[:] pre_theta, pro_theta
    cdef double[:,:] cov
    cdef double pre_lnlike, pro_lnlinke, alpha, uni
    cdef int i

    result = np.zeros((npar+1, nsample))
        
    pre_theta = initheta
    cov = np.eye(npar)*0.1 #np.random.uniform(0,5,npar)
    
    pre_lnlike = lnlike(pre_theta, sigma, dats)
        
    for i in range(nsample):
        pro_theta = np.random.multivariate_normal(pre_theta, cov)
        pro_lnlinke = lnlike(pro_theta, sigma, dats)
        
        alpha = np.min([0, pro_lnlinke -pre_lnlike])
        uni = np.random.uniform()
        if(np.log(uni) < alpha):
            pre_theta = pro_theta
            pre_lnlike = pro_lnlinke
                
        result[0,i] = pre_lnlike
        result[1:,i] = pre_theta[:]
    return result
    '''
    cdef double[:,:] result
    cdef double[:] pre_theta, pro_theta
    cdef double[:,:] cov
    cdef double pre_lnlike, pro_lnlinke, alpha, uni
    cdef int i
    cdef double Ti, Tmax, Tmin, ratio

    Tmax = 100
    Tmin = 1e-5
    ratio = 0.68

    result = np.zeros((npar+1, nsample))
        
    pre_theta = initheta
    cov = np.eye(npar)*0.1 #np.random.uniform(0,5,npar)
    
    pre_lnlike = lnlike(pre_theta, sigma, dats)

    Ti = Tmax
    while(Ti > Tmin):
        for i from 0<= i <nsample:
            pro_theta = np.random.multivariate_normal(pre_theta, cov)
            pro_lnlinke = lnlike(pro_theta, sigma, dats)
            
            if(pro_lnlinke > pre_lnlike):
                pre_theta = pro_theta
                pre_lnlike = pro_lnlinke
            else:
                #alpha = np.min([1., np.exp(pro_lnlinke -pre_lnlike)] )
                alpha = 1.0/(1.0+exp(-(pro_lnlinke -pre_lnlike)/Ti) )
                uni = np.random.uniform()
                if( uni < alpha):
                    pre_theta = pro_theta
                    pre_lnlike = pro_lnlinke
                    
            result[0,i] = pre_lnlike
            result[1:,i] = pre_theta[:]
        Ti = Ti*ratio
    return result

#********************************************************
# pso
#
cpdef double[:,:] pso(lnlike, double[:,:] dats, double[:] sigma, \
        int nsample=500, int npar=2):
    '''
    cpdef double[:,:] pso(lnlike, double[:,:] dats, double[:] sigma, int nsample=500, int npar=2):
    ====================
    cdef double[:,:] Xt, Vt
    cdef double[:,:] Pbest
    cdef double[:] P_fit, Gbest
    cdef double G_fit, tmp
    cdef int npart, i, j, k
    cdef double[:] Xtt, Vtt, xmin, xmax, vmax
    cdef double r1, r2, c1, c2, w, wmax, wmin, t

    npart = 50
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

    xmax[:] = sigma[:]
    xmin[:] = 0.01
    for i from 0<= i < npar:
        vmax[i] = 0.5*(xmax[i] -xmin[i])

    wmax = 0.9
    wmin = 0.1
    t = 0.0
    c1 = 1.193
    c2 = 1.193
    w = wmax


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
