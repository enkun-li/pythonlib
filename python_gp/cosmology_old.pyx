#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name:  
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time:  
#==================================

import numpy as np
import emcee

cdef extern from "math.h":
    double sqrt(double xx)
    double exp(double xx)
    double sinh(double xx)
    double sin(double xx)
    double cos(double xx)
    double fabs(double xx)
    double log10(double xx)
    double log(double xx)

cpdef double lcdm_Hz(double z, double H0=67.27, double omm=0.31, \
        double omk=0.0, double w=-1):
    '''
    lcdm_Hz(double z, double H0=67.27, double omm=0.31, 
            double omk=0.0, double w=-1)
    '''
    cdef double E2z
    E2z = omm*(1+z)**3 + omk*(1+z)**2 +(1-omm-omk)*(1+z)**(3+3*w)
    return H0*sqrt(E2z)

cpdef double omega_mz(double z, double H0=67.27, double omm=0.31, \
        double omk=0.0, double w=-1):
    '''
    omega_mz(double z, double H0=67.27, double omm=0.31,
             double omk=0.0, double w=-1)
    '''
    cdef double h2z, omegam
    h2z = lcdm_Hz(z,H0,omm,omk,w)**2
    omegam = omm*(1+z)**3
    return omegam/h2z

## series to 3-th order
cpdef double Hz_Kin(double z, double[:] theta):
    '''
    Hz_Kin(double z, double[:] theta)
    '''
    cdef double H0,q0,j0,s0
    cdef double y, eps0, eps1, eps2, eps3, hz
    H0,q0,j0,s0 = theta
    y = z/(1.+z)
    eps0 = 1
    eps1 = 1+q0
    eps2 = 1./2.*(2.+j0+2*q0-q0**2)
    eps3 = 1./6.*(6+j0*(3-4*q0)+6*q0 -3*q0**2+3*q0**3-s0)
    #eps4 = -1./24.*(-24 + 4*j0**2- l0 -24*q0 + 12*q0**2 - 12*q0**3 \
    #        + 15*q0**4 + j0* (-12 + 16*q0 - 25*q0**2 ) \
    #        + 4*s0 - 7*q0*s0)
    hz = H0*(eps0 + eps1*y+eps2*y**2+eps3*y**3)
    return hz

cpdef double ComovingDistance_Kin(double z, double[:] theta):
    '''
    ComovingDistance_Kin(double z, double theta)
    '''
    cdef double H0,q0,j0,s0
    cdef double dc
    cdef double y, eps1, eps2, eps3
    cdef double c_unit = 2.99792458e5 # light speed in [km/s]
    H0,q0,j0,s0 = theta
    y = z/(1.+z)
    eps1 = 0.5*(1-q0)
    eps2 = 1./6.*(-j0+3*q0**2-2*q0+2)
    eps3 = 1./24.*(6-j0*(3-10*q0)-6*q0+9*q0**2-15*q0**3+s0)
    #eps4 = 1./120.*(-j0*(105*q0**2-40*q0+12) + 10*j0**2-l0-15*q0*s0 &
    #       +105*q0**4-60*q0**3+36*q0**2-24*q0+4*s0+24)
    
    dc = c_unit*y/H0*(1+eps1*y +eps2*y**2.0+eps3*y**3)
    return dc

cpdef double AngularDiameterDistance_Kin(double z, double[:] theta):
    '''
    AngularDiameterDistance_Kin(double z, double theta)
    '''
    cdef double H0=theta[0], omk=theta[-1]
    cdef double dc, dc_unit, chi, da
    cdef double c_unit = 2.99792458e5 # light speed in [km/s]
    cdef double const_pi = 3.14159265 # pi
    
    dc = ComovingDistance_Kin(z,theta[:4])
    dc_unit = dc/c_unit
    chi = H0*sqrt(fabs(omk))*dc_unit
    
    if(omk>0.0):
        da = c_unit/(1+z)/H0/sqrt(omk)*sinh(chi)
    elif(omk<0.0):
        #if(chi > const_pi):
        #    print('chi >= pi is not allowed')
        #else:
        da = c_unit/(1+z)/H0/sqrt(-omk)*sin(chi)
    else:
        da = dc/(1+z)
    
    return da

cpdef double LuminosityDistance_Kin(double z, double[:] theta):
    '''
    LuminosityDistance_Kin(double z, double theta)
    '''
    cdef double dlum
    dlum = AngularDiameterDistance_Kin(z,theta)
    return dlum

cpdef double BAO_D_V(double z, double[:] theta):
    '''
    BAO_D_V(double z, double theta)
    '''
    cdef double dv, da, hz
    hz = Hz_Kin(z,theta[:4])
    da = AngularDiameterDistance_Kin(z,theta)
    dv = ((1+z)**2*da**2*z/hz)**(1.0/3)
    return dv

## =============================================
## read in CC data
def read_CC_data(filename='CC_31.txt'):
    '''
    read_CC_data(filename='CC_31.txt')
    '''
    cdef double[:,:] dat
    dat = np.loadtxt(filename, unpack=True)
    return dat

cpdef double cc_lnlike(double[:] theta, double[:,:] dats):
    '''
    cc_lnlike(double theta, double[:,:] dats)
    '''
    cdef double[:] zobs, hobs, sigs
    cdef double lnlike, hz_th, chisq, logZero = 1e30
    cdef double difh
    cdef int num
    zobs = dats[0,:]
    hobs = dats[1,:]
    sigs = dats[2,:]
    num = np.shape(zobs)[0]
    lnlike = 0.
    for i from 0 <= i < num:
        hz_th = Hz_Kin(zobs[i],theta)
        if(hz_th < 0):
            lnlike = logZero
            break
        else:
            difh = hz_th - hobs[i]
            lnlike = lnlike + difh**2.0/sigs[i]**2/2.0

    return lnlike

## =============================================
## read in Pantheon data
def read_Pantheon_data(filename='Pantheon.txt'):
    '''
    read_Pantheon_data(filename='Pantheon.txt')
    '''
    cdef double[:,:] dat
    dat = np.loadtxt(filename, unpack=True)
    return dat

def read_Pantheon_covmat(filename='Pantheon_covmat.txt'):
    '''
    read_Pantheon_covmat(filename='Pantheon_covmat.txt')
    '''
    cdef double[:,:]  invcov
    invcov = np.loadtxt(filename, unpack=True)
    return invcov

cpdef double sn_lnlike(double[:] theta, double[:,:] dats, double[:,:] invcov):
    '''
    sn_lnlike(double[:] theta, double[:,:] dats, double[:,:] invcov)
    '''
    cdef double[:] zcmb, zhel, mag, dmag, Acov, diflum
    cdef double da, zi, lumdists, logZero = 1e30
    cdef double amarg_A, amarg_B, amarg_C, lnlike
    cdef int nsn
    zcmb = dats[0,:]
    zhel = dats[1,:]
    mag = dats[2,:]
    dmag = dats[3,:]
    nsn = np.shape(zcmb)[0]
    diflum = np.empty(nsn, dtype=np.float64)
    lnlike = 0

    for i from 0 <= i < nsn:
        zi = zcmb[i]
        da = AngularDiameterDistance_Kin(zi,theta)
        if(da <= 0):
            lnlike = lnlike + logZero
            break
        else:
            lnlike = lnlike + 0
            #lumdists = 5.0*log10((1.0+zcmb)*(1.0+zhel)*da)
            lumdists = 5.0*log10(da)
            diflum[i] = lumdists - mag[i]

    if(lnlike >= 1e30):
        lnlike = logZero
    else:
        Acov = np.matmul(diflum, invcov)
        amarg_A = np.dot(Acov, diflum)
        amarg_B = np.sum(Acov)
        amarg_C = np.sum(invcov)
        lnlike = amarg_A - amarg_B**2/amarg_C
    return lnlike

cpdef double sn_mag(double[:] theta, double[:,:] dats, double[:,:] invcov):
    '''
    sn_mag(double theta, double[:,:] dats, double[:,:] invcov)
    '''
    cdef double[:] zcmb, zhel, mag, dmag, Acov, diflum
    cdef double da, zi, lumdists
    cdef double amarg_B, amarg_C, lnlike
    cdef int nsn
    zcmb = dats[0,:]
    zhel = dats[1,:]
    mag = dats[2,:]
    dmag = dats[3,:]
    nsn = np.shape(zcmb)[0]
    diflum = np.empty(nsn, dtype=np.float64)

    for i from 0 <= i < nsn:
        zi = zcmb[i]
        da = AngularDiameterDistance_Kin(zi,theta)
        lumdists = 5.0*log10(da)
        diflum[i] = lumdists - mag[i]

    Acov = np.matmul(diflum, invcov)
    amarg_B = np.sum(Acov)
    amarg_C = np.sum(invcov)
    return amarg_B/amarg_C

## ==========================================
# priors on parameters
cpdef double cc_lnprior(double[:] theta):
    '''
    cc_lnprior(double[:] theta)
    '''
    cdef double H0,q0,j0,s0
    H0,q0,j0,s0 = theta
    if( 50 < H0 < 80 and -1.0 < q0 < 0.0 and \
            -20 < j0 < 20 and -200 < s0 < 100 \
            ):
        return 0.0
    return -1e30

# probability
cpdef double cc_lnprob(double[:] theta, double[:] zob, double[:] Hob, \
        double[:] sigob):
    '''
    cc_lnprob(double[:] theta, double[:,:] cc_dat)
    '''
    cdef double[:,:] cc_dat
    cdef double lp = cc_lnprior(theta)    
    cc_dat = np.array([zob, Hob, sigob])
    if(lp <= -1e30):
        return -1e30
    return lp + cc_lnlike(theta,cc_dat)

# priors on parameters
cpdef double sn_lnprior(double[:] theta):
    '''
    sn_lnprior(double H0, double q0, double j0,
                double s0, double omk)
    '''
    cdef double H0,q0,j0,s0,omk
    H0,q0,j0,s0,omk = theta
    if(50 < H0 < 80 and -1.0 < q0 < 0.0 and \
            -20 < j0 < 20 and -200 < s0 < 100 \
            and -1 < omk < 1):
        return 0.0
    return -1e30

# probability
cpdef double sn_lnprob(double[:] theta, double[:,:] sn_dat, double[:,:] sn_cov):
    '''
    sn_lnprob(double theta, double[:,:] sn_dat, double[:,:] sn_cov)
    '''
    cdef double lp = sn_lnprior(theta)
    if(lp <= -1e30):
        return -1e30
    return lp + sn_lnlike(theta, sn_dat, sn_cov)

## ==========================================
# run a mcmc
def cc_runmcmc(int ndim, int nwalkers, int nsteps, double[:,:] cc_dat,\
        double[:] theta_c):
    '''
    cc_runmcmc(int ndim, int nwalkers, int nsteps, double cc_dat, double theta_c)
    '''
    cdef double[:,:] pos
    
    pos = theta_c + 1e-3*np.random.randn(nwalkers, ndim)
    sampler_cc = emcee.EnsembleSampler(nwalkers,ndim,cc_lnprob, args=cc_dat)
    sampler_cc.run_mcmc(pos, nsteps, progress=True)
    return sampler_cc

## ==========================================
def array_func(func, theta, array):
    '''
    output array of function theta in the array.
    '''
    return np.array([func(xx, theta) for xx in array])

## ============================================
## random choice of priors for parameters
def randomchoice():
    '''
    output 5 random double parameters under uniform distributions.
    '''
    cdef double H0, q0, j0, s0, omk
    H0 = np.random.uniform(50,80)
    q0 = np.random.uniform(-1,0)
    j0 = np.random.uniform(-20,20)
    s0 = np.random.uniform(-200,100)
    omk = np.random.uniform(-1,1)
    return np.array([H0,q0,j0,s0,omk])

## ============================================================================
def findpar(int N):
    cdef double[:] theta
    cdef double hz, da
    cdef int i, j
    cdef double[:,:] thetas
    
    thetas = np.empty((5,N), dtype=np.float64)
    
    for i from 0 <= i < N:
        theta = randomchoice()
        hz = Hz_Kin(2.0, theta[:4])
        da = AngularDiameterDistance_Kin(2.0, theta)
        if(hz > 0 and da > 0):
            for j from 0 <= j < 5:
                thetas[j,i] = theta[j]
        else:
            for j from 0 <= j < 5:
                thetas[:,i] = 0.0
    return thetas

