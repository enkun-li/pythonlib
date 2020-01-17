#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: rec_cosmology.pyx
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2018-12-25 21:29:48
#==================================

cimport cython

import numpy as np
import scipy.optimize as opt
from scipy import integrate
from gapp import dgp

from gaussian_process import recfuncpar
from gaussian_process import gaussian_process as gps

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

#******************************************
# reconstract cosmology background
#
cdef class rec_cosmology:
    cdef double omk, H0
    cdef double c_light, curv2, curv
    cdef int Ksign
    cdef char* flat
    
    def __cinit__(self, double omk, double H0):
        '''
        reconstract cosmology background:
        ==================================
        input: omega_K & H0
        '''
        self.omk = omk
        self.H0 = H0
        self.c_light = 2.99792458e8
        self.curv2 = omk*H0**2/(self.c_light)**2
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
        print('This is a %s universe.'%self.flat)
        
    def setup_rec_cosmo(self, filename):
        '''
        setup_rec_cosmo(filename)
        ==============================
        return sigf, lenf
        '''
        cdef double[:] initheta, theta
        
        initheta = np.array([100,2], dtype=np.float64)
        pars = recfuncpar(filename, initheta)
        theta = pars.theta()
        
        gpsrec = gps(filename)
        gpsrec.setup_gp(theta[0], theta[1])
        return gpsrec

    
    
    

