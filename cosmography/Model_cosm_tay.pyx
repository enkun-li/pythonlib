#!/home/ekli/.local/anaconda3/bin/python
#-*- coding: utf-8 -*-  
#==================================
# File Name: Model_cosm_tay.pyx
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-10-22 15:13:39
#==================================

cimport cython
import numpy as np
from scipy.integrate import quad
import sys

cdef extern from "math.h":
    double pow(double x, double y)
    double sqrt(double x)
    double sinh(double x)
    double sin(double x)
    double fabs(double x)

cdef class Model_cosm_tay:
    cdef double H0, q0, j0, s0, l0, omk
    cdef double fy1, fy2, fy3, fy4
    cdef double omksq
    cdef double c
    cdef int case #[i: Ti]

    def __cinit__(self, double H0=70.0, double q0=-0.55, double j0=1.0, \
            double s0=-0.35, double l0=3.115, double omk=0, int case=3):
        self.H0 = H0
        self.q0 = q0
        self.j0 = j0
        self.s0 = s0
        self.l0 = l0
        self.omk = omk
        self.omksq = omk*omk
        self.case = case
        self.c = 2.99792458e5 #[km/s]
        #===================
        self.fy1 = (1-q0)
        self.fy2 = (2-j0 -2*q0 +3*q0**2 +omk)
        self.fy3 = (6-j0 * (3-10 * q0)-6*q0 +9 *q0**2 -15 *q0**3 \
                +s0+6 *(1-q0) *omk)
        self.fy4 = (24+10 *j0**2 -l0-24 *q0 +36 *q0**2 -60*q0**3 \
                +105*q0**4-j0 *(12-40 *q0 +105* q0**2)+4*s0-15*q0*s0 \
                -5*(-7+2*j0+10*q0-9*q0**2) *omk +self.omksq)
        return

    ## ===================================    
    cpdef double taylor_cg_dy2(self, double z):
        cdef double y = z/(1+z)
        return y + (self.fy1 * y**2)/2

    cpdef double taylor_cg_dy3(self, double z):
        cdef double y = z/(1+z)
        return y + (self.fy1 * y**2)/2 + (self.fy2 * y**3)/6
    
    cpdef double taylor_cg_dy4(self, double z):
        cdef double y = z/(1+z)
        return y + (self.fy1 * y**2)/2 + (self.fy2 * y**3)/6 \
                + (self.fy3 * y**4)/24
                
    cpdef double taylor_cg_dy5(self, double z):
        cdef double y = z/(1+z)
        return y + (self.fy1 * y**2)/2 + (self.fy2 * y**3)/6 \
                + (self.fy3 * y**4)/24 + (self.fy4 * y**5)/120

    ## -----------------------------------
    cpdef double taylor_cg_pd_dy2(self, double z):
        cdef double y = z/(1+z)
        return 1 + (self.fy1 * y)
    
    cpdef double taylor_cg_pd_dy3(self, double z):
        cdef double y = z/(1+z)
        return 1 + (self.fy1 * y) + (self.fy2 * y**2)/2
    
    cpdef double taylor_cg_pd_dy4(self, double z):
        cdef double y = z/(1+z)
        return 1 + (self.fy1 * y) + (self.fy2 * y**2)/2 \
                + (self.fy3 * y**3)/6

    cpdef double taylor_cg_pd_dy5(self, double z):
        cdef double y = z/(1+z)
        return 1 + (self.fy1 * y) + (self.fy2 * y**2)/2 \
                + (self.fy3 * y**3)/6 + (self.fy4 * y**4)/24
    
    ## ===================================
    cpdef double Hofz(self, double z):        
        cdef double hz
        hz = 0
        if(self.case == 2):
            hz = (1+z)**2 * self.H0/self.taylor_cg_pd_dy2(z) \
                    * sqrt(1+self.omk*(self.taylor_cg_dy2(z))**2)
        elif(self.case == 3):
            hz = (1+z)**2 * self.H0/self.taylor_cg_pd_dy3(z) \
                    * sqrt(1+self.omk*(self.taylor_cg_dy3(z))**2)
        elif(self.case == 4):
            hz = (1+z)**2 * self.H0/self.taylor_cg_pd_dy4(z) \
                    * sqrt(1+self.omk*(self.taylor_cg_dy4(z))**2)
        elif(self.case == 5):
            hz = (1+z)**2 * self.H0/self.taylor_cg_pd_dy5(z) \
                    * sqrt(1+self.omk*(self.taylor_cg_dy5(z))**2)
        else:
            sys.exit("No such model %s! [Model_cosm_tay]"%self.case)
        return hz

    cpdef double __dz(self, double z):
        cdef double dz
        dz = 0
        if(self.case == 2):
            dz = self.taylor_cg_dy2(z)
        elif(self.case == 3):
            dz = self.taylor_cg_dy3(z)
        elif(self.case == 4):
            dz = self.taylor_cg_dy4(z)
        elif(self.case == 5):
            dz = self.taylor_cg_dy5(z)
        else:
            sys.exit("No such model %s! [Model_cosm_tay]"%self.case)
        return dz

    cpdef double distance_DM(self, double z):
        return self.c/self.H0 * self.__dz(z)

    cpdef double distance_DA(self, double z):
        return 1/(1+z)*self.c/self.H0 * self.__dz(z)

    cpdef double distance_DL(self, double z):
        return (1+z) * self.c/self.H0 * self.__dz(z)

    cpdef double distance_DV(self, double z):
        cdef double dv3
        dv3 = (self.distance_DM(z))**2 * self.c*z/self.Hofz(z)
        return (dv3)**(1.0/3.0)

