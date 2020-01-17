#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: likelihood.pyx
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time:  
#==================================

cimport cython

import numpy as np

cdef class sn_likelihood:
    #cdef double[:,:] sndata, rec_dc
    cdef double[:,:] syscov, rec_covdc, invcov
    cdef double[:] zcmb, snm, sndm, recdc
    cdef int nsn
    cdef double H0, omk, fk, c_light
    cdef double flatK

    def __cinit__(self,datafile,sysfile,recdcfile,reccovfile):
        '''
        initial sn_likelihood
        =======================================================
        datafile: the filename of SN data
        sysfile:  the filename of system covariance
        recdcfile: the filename of reconstract d_c at z(SN)
        reccovfile: the filename of covariance of cov(dc(zi),dc(zj))
        '''
        cdef double[:,:] sndata, rec_dc

        self.c_light = 2.99792458e5 # [km/s]
        self.flatK = 1e-8

        sndata = np.load(datafile)
        self.syscov = np.load(sysfile)
        rec_dc = np.load(recdcfile)
        self.rec_covdc = np.load(reccovfile)

        self.zcmb = sndata[0]
        self.snm = sndata[1]
        self.sndm = sndata[2]
        self.recdc = rec_dc[1]

        self.nsn = np.shape(sndata)[1]
        self.invcov = np.zeros((self.nsn, self.nsn), dtype=np.float64)
        #print('='*50)
        #print('Initial sn likelihood')
        #print('The number of SN is: %s'%self.nsn)

    def setup_pars(self, double H0, double omk):
        '''
        setup the model parameters
        ================================
        H0: hubble parameter
        omk: spatial curvature
        '''
        self.H0 = H0
        self.omk = omk
        if(np.abs(omk) < self.flatK ):
            self.fk = 1.0
        else:
            self.fk = H0*np.sqrt(np.abs(omk))/self.c_light

        self.invcov = self.InvMcov()
        #print('H0 is: %s'%self.H0)
        #print('Omk is: %s'%self.omk)
        #print('fk is: %s'%self.fk)
        
    cpdef double __sinn(self, double chi):
        '''
        sinn(chi):
        =========================
        sinh(chi) for omk > 0
        chi       for omk = 0
        sin(chi)  for omk < 0
        '''
        cdef double sinnchi
        if(self.omk > self.flatK):
            sinnchi = np.sinh(chi)
        elif(self.omk < -self.flatK):
            sinnchi = np.sin(chi)
        else:
            sinnchi = chi
        return sinnchi
    
    cpdef double __cosn(self, double chi):
        '''
        cosn(chi):
        =========================
        cosh(chi) for omk > 0
        chi       for omk = 0
        cos(chi)  for omk < 0
        '''
        cdef double cosnchi
        if(self.omk > self.flatK):
            cosnchi = np.cosh(chi)
        elif(self.omk < -self.flatK):
            cosnchi = np.cos(chi)
        else:
            cosnchi = chi/chi
        return cosnchi

    cpdef double[:,:] InvMcov(self):
        '''
        InvMcov
        ==================
        '''
        cdef double[:,:] covuu, covtot, invcov
        cdef double[:] chi, vec1, vec2, vec
        cdef int i, j
        
        covtot = np.zeros((self.nsn,self.nsn))
        invcov = np.zeros((self.nsn, self.nsn))
        covuu = np.zeros((self.nsn, self.nsn), dtype=np.float64)
        chi = np.zeros(self.nsn, dtype=np.float64)
        vec = np.zeros(self.nsn, dtype=np.float64)
        
        for i from 0 <= i < self.nsn:
            for j from 0 <= j < self.nsn:
                if(i==0):
                    chi[j] = self.fk*self.recdc[j]
                    vec[j] = 5.0/np.log(10.0)*self.fk*self.__cosn(chi[j])/self.__sinn(chi[j])
                covuu[i,j] = vec[i]*vec[j]*self.rec_covdc[i,j]
                covtot[i,j] = covuu[i,j] + self.syscov[i,j]
            covtot[i,i] = covtot[i,i] + self.sndm[i]

        invcov = np.linalg.inv(covtot)

        return invcov

    
    cpdef double[:] amag(self, double H0, double omk):
        '''
        lnlike_sn
        =======================
        H0:
        omk:
        '''
        cdef double[:] difmu, vec, amag_M
        cdef double A, B, C, chisq, dlu
        cdef int i

        self.setup_pars(H0,omk)

        difmu = np.zeros(self.nsn)
        vec = np.zeros(self.nsn)
        amag_M = np.zeros(3)

        for i from 0 <= i < self.nsn:
            dlu = (1+self.zcmb[i])/self.fk*self.__sinn(self.fk*self.recdc[i])
            difmu[i] = self.snm[i] -5.0/np.log(10)*np.log(dlu)
        vec = np.matmul(difmu,self.invcov)
        A = np.dot(vec,difmu)
        B = np.sum(vec)
        C = np.sum(self.invcov)
        amag_M[0] = A
        amag_M[1] = B
        amag_M[2] = C
        return amag_M 

    cpdef double lnlike_sn(self, double H0, double omk):
        '''
        lnlike_sn
        =======================
        H0:
        omk:
        '''
        cdef double[:] mag
        cdef double chisq
        cdef int i

        mag = self.amag(H0,omk)

        chisq = mag[0] -mag[1]**2/mag[2] +np.log(mag[2]/2.0/np.pi)
        chisq = -chisq/2.0
        return chisq


    cpdef double lnprior(self, double[:] theta):
        '''
        lnprior(theta)
        ===========================
        '''
        cdef double H0, omk
        H0, omk = theta
        if(60 < H0 < 80 and -1<omk<1):
            return 0.0
        return -np.inf

    cpdef double lnprob(self,double[:] theta):
        '''
        lnprob(theta)
        =========================
        '''
        cdef double H0, omk, lp
        lp = self.lnprior(theta)
        H0,omk = theta
        if(not np.isfinite(lp)):
            return -np.inf
        return lp+self.lnlike_sn(H0, omk)
