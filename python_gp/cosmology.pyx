#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name:  
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time:  
#==================================

cimport cython
import numpy as np
import scipy.optimize as opt
from scipy import integrate

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

#*****************************
# cosmology background
#
cdef class cosmology_bk:
    cdef double omm, omk, w, H0, omx
    cdef double c_light, curv2, curv
    cdef int Ksign
    cdef char* flat

    def __cinit__(self, double omm, double omk, double w, double h0):
        self.omm = omm
        self.omk = omk
        self.omx = 1.0-omm-omk
        self.w = w
        self.H0 = h0
        self.c_light = 2.99792458e5 #[km/s]
        self.curv2 = omk*h0**2/(self.c_light)**2
        if(self.curv2 > 1e-8):
            # open
            self.curv = sqrt(self.curv2)
            self.Ksign = -1
            self.flat = 'open'
        elif(self.curv2 < -1e-8):
            # closed
            self.curv = sqrt(-self.curv2)
            self.Ksign = 1
            self.flat = 'closed'
        else:
            # flat
            self.curv = 1.0
            self.Ksign = 0
            self.flat = 'flat'
        #print('This is a %s universe.'%self.flat)

    
    cpdef double Hofz(self, double z):
        '''
        Hofz(z)
        ==============
        in [Mpc]
        '''
        cdef double hz
        hz = self.omm*(1.0+z)**3.0 +self.omx*(1+z)**(3.0*(1+self.w)) \
                +self.omk*(1+z)**2.0
        return sqrt(hz)*self.H0

    cpdef double sigmaHz(self, double z, double sigh0, double sigomm, \
            double sigomk):
        '''
        sigmaHz(z,sigh0,sigomm,sigomk)
        =====================================
        sig^2 = \sum \pd H/\pd p_i \pd H/\pd p_j cov(p_i,p_j)
        '''
        cdef double sig, sig1, sig2, sig3
        sig1 = self.Hofz(z)**2/self.H0**2*sigh0**2
        sig2 = self.H0**4/4.0/self.Hofz(z)**2*( ((1+z)**3-1)**2*sigomm**2)
        sig3 = self.H0**4/4.0/self.Hofz(z)**2*( ((1+z)**2-1)**2*sigomk**2)
        sig = sqrt(sig1+sig2+sig3)
        return sig

    cpdef double __dtauda(self, double z):
        '''
        dtauda(z)
        ==================
        1/H(z)
        '''
        return 1.0/self.Hofz(z)
    
    cpdef double omegam(self, double z):
        '''
        omegam(z)
        ==================
        '''
        cdef double ff
        ff = self.omm*(1+z)**3.0/self.Hofz(z)
        return ff

    cpdef double omegax(self, double z):
        '''
        omegax(z)
        ===============
        '''
        cdef double ff
        ff = self.omx*(1+z)**(3*(1+self.w))/self.Hofz(z)
        return ff

    cpdef double comovingdistance(self, double z):
        '''
        comovingdistance(z)
        ========================
        dc in Mpc
        '''
        cdef double dc, err
        if(z==0.0):
            dc = 0.0
        else:
            dc, err = integrate.quad(self.__dtauda, 0.0, z)
        return self.c_light*dc

    cpdef double __dtaudasqsig(self, double z, double sigh0, double sigomm, double sigomk):
        '''
        __dtaudasqsig(z,sigh0,sigomm,sigomk)
        =====================================
        1/dtauda(z)
        '''
        cdef double fsh, ff
        fsh = self.sigmaHz(z,sigh0,sigomm,sigomk)
        ff = self.__dtauda(z)**2.0*fsh
        return ff

    cpdef double sigma_dc(self, double z, double sigh0, double sigomm, double sigomk):
        '''
        sigma_dc(z,sigh0,sigomm,sigomk)
        ====================================
        sigma_dc = c \int_0^z sigmaHz(x)/H(x)^2 dx
        '''
        cdef double ddci, err, sigh, sigdc
        if(z==0.0):
            sigdc = 0.0
        else:
            ddci, err = integrate.quad(self.__dtaudasqsig,0,z, args=(sigh0,sigomm,sigomk))
            sigdc = fabs(ddci)*self.c_light
        return sigdc

    cpdef double __sinn(self, double chi):
        '''
        __sinn(chi)
        ==========================
        sinh(chi) for open
        sin(chi)  for closed
        chi       for flat
        '''
        
        if(self.Ksign == -1):
            return sinh(chi)
        elif(self.Ksign == 1):
            return sin(chi)
        else:
            return chi
    
    cpdef double __cosn(self, double chi):
        '''
        __cosn(chi)
        ==========================
        cosh(chi) for open
        cos(chi)  for closed
        1         for flat
        '''
        
        if(self.Ksign == -1):
            return cosh(chi)
        elif(self.Ksign == 1):
            return cos(chi)
        else:
            return 1.0

    cpdef double comovingangulardistance(self, double z):
        '''
        comovingangulardistance(z)
        ================================
        d_M = (1+z)dA in [Mpc]
        '''
        cdef double dM, chi
        chi = self.curv*self.comovingdistance(z)
        dM = 1.0/self.curv*self.__sinn(chi)
        return dM
   
    cpdef double sigma_fc(self,double sigh0, double sigomk):
        '''
        sigma_fc(sigh0,sigomk)
        ================================
        fc = H0\sqrt{omk}/c
        sig_fc2 = (\sqrt(omk)/c)^2*sigh0^2 
            + (H0/c/2/\sqrt{omk})^2 *sigomk^2
        '''
        cdef double sigfc2, h0, omk, c, sigfc
        h0 = self.H0
        omk = fabs(self.omk)
        c = self.c_light
        if(omk < 1e-8):
            sigfc = 0.0
        else:
            sigfc2 = (sqrt(omk)/c)**2* sigh0**2.0 +(h0/c/2.0/sqrt(omk) )**2*sigomk**2
            sigfc = sqrt(sigfc2)
        return sigfc


    cpdef double sigma_dm(self, double z, double sigh0, double sigomm, double sigomk):
        '''
        sigma_dm(z,sigh0,sigomm,sigomk)
        ================================
        sigma_dm = cosn(fc dc) sigma_dc
        '''
        cdef double dc, ddM, chi, dmdc, dmdf
        dc = self.comovingdistance(z)
        chi = self.curv * dc
        dmdc = self.__cosn(chi)*self.sigma_dc(z,sigh0,sigomm,sigomk)
        dmdf = (dc/self.curv * self.__cosn(chi) -1.0/self.curv**2 *self.__sinn(chi)) \
                * self.sigma_fc(sigh0,sigomk)
        ddM = sqrt(dmdc**2 +dmdf**2)
        return ddM


    cpdef double angulardistance(self, double z):
        '''
        angulardistance(z)
        ============================
        d_A = d_M/(1+z)
        '''
        return self.comovingangulardistance(z)/(1.0+z)

    cpdef double luminositydistance(self, double z):
        '''
        luminositydistance(z)
        =============================
        d_L = (1+z) d_M = (1+z)^2 d_A
        '''
        return (1+z)*self.comovingangulardistance(z)

    cpdef double distancemodulus(self, double z, double mu0=25.0):
        '''
        distancemodulus(z)
        ==============================
        std: mu = 5\log_{10}[dl/(10pc)] = 5 lg[dl/Mpc] +25
        in [Mpc]
        '''
        cdef double mu, dl
        dl = self.luminositydistance(z)
        mu = 5.0*log10(dl)+mu0
        return mu

    cpdef double sigma_muth(self, double z, double sigh0,double sigomm,double sigomk):
        '''
        sigma_muth(z,sigh0,sigomm,sigomk)
        =====================================
        sigma_muth = 5/ln10/dm*sigma_dm
        '''
        cdef double dm, sigdm,sigmu
        if(z==0):
            sigmu = 0.0
        else:
            dm = self.comovingangulardistance(z)
            sigdm = self.sigma_dm(z,sigh0,sigomm,sigomk)
            sigmu = 5.0/log(10.0)/dm*sigdm
        return sigmu

