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
    cdef double[:,:] sndata, rec_dc
    cdef double[:,:] syscov, rec_covdc
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
        self.c_light = 2.99792458e5 # [km/s]
        self.flatK = 1e-8

        self.sndata = np.load(datafile)
        self.syscov = np.load(sysfile)
        self.rec_dc = np.load(recdcfile)
        self.rec_covdc = np.load(reccovfile)
        self.nsn = np.shape(self.sndata)[1]
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
        #print('H0 is: %s'%self.H0)
        #print('Omk is: %s'%self.omk)
        #print('fk is: %s'%self.fk)

    cpdef double lnlike_sn(self, double H0, double omk):
        '''
        lnlike_sn
        =======================
        H0:
        omk:
        '''
        cdef double[:,:] invcov
        cdef double[:] difmu, vec
        cdef double A, B, C, chisq
        cdef int i

        self.setup_pars(H0, omk)

        difmu = np.zeros(self.nsn)
        vec = np.zeros(self.nsn)
        invcov = self.__InvMcov()

        for i from 0 <= i < self.nsn:
            difmu[i] = self.sndata[1,i] -5.0/np.log(10)*np.log(self.rec_dc[1,i])
        vec = np.matmul(difmu,invcov)
        A = np.dot(vec,difmu)
        B = np.sum(vec)
        C = np.sum(invcov)
        chisq = A -B**2/C        
        return chisq
    
    cpdef double[:] amag(self, double H0, double omk):
        '''
        lnlike_sn
        =======================
        H0:
        omk:
        '''
        cdef double[:,:] invcov
        cdef double[:] difmu, vec, amag_M
        cdef double A, B, C, chisq
        cdef int i

        self.setup_pars(H0,omk)

        difmu = np.zeros(self.nsn)
        vec = np.zeros(self.nsn)
        amag_M = np.zeros(3)
        invcov = self.__InvMcov()        

        #print('calculate the diff between th and obs')

        for i from 0 <= i < self.nsn:
            difmu[i] = self.sndata[1,i] -5.0/np.log(10)*np.log(self.rec_dc[1,i])
        vec = np.matmul(difmu,invcov)
        A = np.dot(vec,difmu)
        B = np.sum(vec)
        C = np.sum(invcov)
        chisq = A -B**2/C 
        amag_M[0] = A
        amag_M[1] = B
        amag_M[2] = C
        return amag_M

    cpdef double[:,:] __InvMcov(self):
        '''
        InvMcov
        ==================
        '''
        cdef double[:,:] covuu, covtot, invcov
        cdef double[:] chi, vec1, vec2, vec
        cdef int i, j
        
        covuu = np.zeros((self.nsn, self.nsn), dtype=np.float64)
        covtot = np.zeros((self.nsn, self.nsn), dtype=np.float64)
        chi = np.zeros(self.nsn, dtype=np.float64)
        #vec1 = np.zeros(self.nsn, dtype=np.float64)
        #vec2 = np.zeros(self.nsn, dtype=np.float64)
        vec = np.zeros(self.nsn, dtype=np.float64)
        
        for i from 0 <= i < self.nsn:
            #chi[i] = self.fk*self.rec_dc[1,i]
            #vec1[i] = 5.0/np.log(10.0)*self.fk*self.__cosn(chi[i])/self.__sinn(chi[i])
            for j from 0 <= j < self.nsn:
                if(i==0):
                    chi[j] = self.fk*self.rec_dc[1,j]
                    vec[j] = 5.0/np.log(10.0)*self.fk*self.__cosn(chi[j])/self.__sinn(chi[j])
                covuu[i,j] = vec[i]*vec[j]*self.rec_covdc[i,j]
                covtot[i,j] = covuu[i,j] + self.syscov[i,j]
            covtot[i,i] = covtot[i,i] + self.sndata[2,i]

        invcov = np.linalg.inv(covtot)

        return invcov

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
            sinnchi = np.sinh(self.fk*chi)
        elif(self.omk < -self.flatK):
            sinnchi = np.sin(self.fk*chi)
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
            cosnchi = np.cosh(self.fk*chi)
        elif(self.omk < -self.flatK):
            cosnchi = np.cos(self.fk*chi)
        else:
            cosnchi = 1.0 
        return cosnchi

