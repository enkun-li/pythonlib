#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: equations.pyx
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-04-14 07:38:28
#==================================

import sys
if(not sys.path[0][-10:-1] == 'pythonlib'):
    sys.path.insert(0, '/home/ekli/myworks/pythonlib/')
from constants import constants as const
import numpy as np

cdef class LambdaGeneral:
    '''
    LambdaGeneral()
    ==============================
    General cosmological model, for dark matter and dark energy, etc.
    param H0: Hubble constant
    param omm: present dimensionless energy density of matter (include dark matter,
        baryon)
    param omr: present dimensionless energy density of radation
    param omk: present dimensionless energy density of spatial curvature
    param omx: present dimensionless energy density of dark energy
    param w_lam: p/rho for the dark energy (assumed constant -1)
    param wa_ppf: for wCDM model
    param cs2_lam: sound speed of dark energy
    param nps(:): array like, model parameters
    param mod: different dark energy models or 
        0: cosmological constant model (LCDM)
        1: constant EOS model (wCDM)
        2: CPL model (wwaCDM)
        3: early dark energy model
    '''
    cdef double w_lam, wa_ppf, cs2_lam
    cdef double H0, omx, omm, omr, omk#, omr
    cdef double[:] nps
    cdef int mod
    
    def __cinit__(self, double H0=67.27, double omm=0.31, double omr=0.0,
            double omk=0.0, double omx=0.69, nps=None,
            double w_lam=-1, double wa_ppf=0, double cs2_lam=1,
            int mod=0):
        self.H0 = H0
        self.w_lam = w_lam
        self.wa_ppf = wa_ppf
        self.cs2_lam = cs2_lam
        self.omm = omm
        self.omk = omk
        self.omx = omx
        if(nps is not None):
            self.nps = nps
        else:
            self.nps = np.zeros(10)
        self.mod = mod

    cpdef double[:] initial_ede(self, double a):
        ''' 
        initial_ede(a)
        ===============================
        Initial early dark energy model.
        param a: scale factor
        return omxi: array like parameters for ede model parameters
        '''
        cdef double w0, omxe, omm0, omx0
        cdef double[:] omxi
        omxi = np.zeros(10)
        w0 = self.nps[0]
        omxe = self.nps[1]
        omm0 = self.omm
        omx0 = self.omm*self.omx/(1-self.omx)
        
        omxi[1] = omxe*(1 -a**(-3*w0))
        omxi[2] = omx0 + omm0*a**(3*w0) #+ omr0*a**(3*w0-1)
        omxi[3] = omx0 - omxi[1]
        
        omxi[4] = 3*w0*omxe*a**(-3*w0-1)
        omxi[5] = 3*w0*omm0*a**(3*w0-1) #+(3*w0-1)*omr0*a**(3*w0-2)
        omxi[6] = -omxi[4]

        omxi[7] = -3*w0*(3*w0 +1)*omxe*a**(-3*w0-2)
        omxi[8] = 3*w0*(3*w0-1)*omm0*a**(3*w0-2) # \
        #    +(3*w0-1)*(3*w0-2)*omr0*a**(3*w0-3)i
        omxi[9] = -omxi[7]
        return omxi
        
    cpdef double omega_x(self, double a):
        '''
        omega_x(a)
        ============================
        Dimensionless energy density of ede
        omx = f3/f2+f1
        param a: scale factor
        '''
        cdef double[:] omxi
        cdef omx
        omxi = self.initial_ede(a)
        omx = omxi[3]/omxi[2] +omxi[1]
        return omx

    cpdef double domxda(self, double a): #domx/da
        '''
        domxda(a)
        =================================
        domx/da
        param a: scale factor
        '''
        cdef double[:] omxi
        cdef double domx
        
        omxi = self.initial_ede(a)
        domx = omxi[6]/omxi[2] - omxi[3]*omxi[5]/omxi[2]**2 +omxi[4]
        return domx

    cpdef double d2omxda2(self, double a): #d2omx/da2
        '''
        d2omxda2(a)
        =============================
        d2omx/da2
        param a: scale factor
        '''
        cdef double[:] omxi
        cdef double d2omx

        omxi = np.zeros(10)
        omxi = self.initial_ede(a)
        d2omx = omxi[9]/omxi[2] -2*omxi[6]*omxi[5]/omxi[2]**2 \
                - omxi[3]*omxi[8]/omxi[2]**2+2*omxi[3]*omxi[5]**2/omxi[2]**3 \
                +omxi[7]
        return d2omx

    cpdef double E2ofa_ede(self, double a):
        '''
        E2ofa_ede(a)
        ===========================
        E2=H^2/H0^2 = (omm a+ omr)a^{-4}/(1-Omeaga_x(a))
        param a: scale factor
        '''
        cdef double E2
        E2 = self.omk*a**2 +self.omm*a+self.omr
        E2 = E2/(1-self.omega_x(a))/a**4
        return E2

    cpdef double E2a4_node(self, a):
        '''
        E2a4_node(a)
        ========================
        without dark energy * a^4
        '''
        cdef double E2
        E2 = self.omk*a**2 +self.omm*a+self.omr
        return E2

    cpdef double Pa4_node(self, a):
        '''
        Pa4_node(a)
        ========================
        pressure without dark energy * a^4
        '''
        cdef double P
        P = -1.0/3.0*self.omk*a**2 +1.0/3.0*self.omr
        return P
    
    cpdef double w_de(self, a):
        '''
        w_de(a)
        ================================
        Equation of state parameter of dark energy
        '''
        cdef double w

        if(self.mod == 0):
            w = self.w_lam
        elif(self.mod == 1):
            w = self.w_lam
        elif(self.mod == 2):
            w = self.w_lam + self.wa_ppf*(1-a)
        elif(self.mod == 3):
            w = self.Pa4_node(a)/self.E2a4_node(a) \
                    - 1.0/3.0*self.domxda(a)/self.omega_x(a)/(1-self.omega_x(a))
        else:
            print("No such model, check it in equations.pyx [w_de]")
            sys.exit(0)

    cpdef double dtauda(self, a):
        '''
        dtauda(a)
        =========================
        dtauda = c/H(a)/a^2 [light speed c in km/s]
        param a: scale factor
        '''
        cdef double E2a4, dtda
        E2a4 = self.omk*a**2 +self.omm*a +self.omr

        if(self.mod == 0):
            E2a4 += self.omx*a**4
        elif(self.mod == 1):
            E2a4 += self.omx*a**(1-3*self.w_lam)
        elif(self.mod == 2):
            E2a4 += self.omx*a**(1-3*(self.w_lam +self.wa_ppf)) \
                    *np.exp(-3*self.wa_ppf*(1-a))
        elif(self.mod == 3):
            E2a4 = E2a4/(1-self.omega_x(a))
        else:
            print("No such model, check it in equations.pyx [dtauda]")
            sys.exit(0)

        dtda = const['c']/1000.0/self.H0/np.sqrt(E2a4)
        return dtda


