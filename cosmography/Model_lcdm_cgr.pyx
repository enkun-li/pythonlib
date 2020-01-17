#!/home/ekli/.local/anaconda3/bin/python
#-*- coding: utf-8 -*-  
#==================================
# File Name: Model_lcdm_cgr.pyx
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-10-22 14:59:12
#==================================

cimport cython
import numpy as np
from scipy.integrate import quad

cdef extern from "math.h":
    double pow(double x, double y)
    double sqrt(double x)
    double sinh(double x)
    double sin(double x)
    double fabs(double x)

cdef class Model_lcdm_cgr:
    cdef double H0, q0, j0, s0, l0, omk
    cdef double c

    def __cinit__(self, double H0, double q0, double j0, double s0, double l0, double omk=0):
        self.H0 = H0
        self.q0 = q0
        self.j0 = j0
        self.s0 = s0
        self.l0 = l0
        self.omk = omk
        self.c = 2.99792458e5 #[km/s]
        return

    cpdef double Hofz(self, double z):
        cdef double hz2
        cdef double f1, f2, f3, f4
        f1 = 2*(1+self.q0)
        f2 = 1+self.j0+2*self.q0
        f3 = -1/3.0*(self.j0*self.q0 +self.s0)
        f4 = 1/12.0*(self.l0 -self.j0*self.j0 \
                +(4-3*self.q0)*(self.q0*self.j0 +self.s0))
        
        return self.H0*sqrt(1 +f1*z +f2*pow(z,2) +f3*pow(z,3) +f4*pow(z,4) )

    cpdef double dtauda(self, double z):
        return self.H0/self.Hofz(z)

    cpdef double __dc(self, double z):
        cdef double dc, err
        dc, err = quad(self.dtauda, 0.0, z)
        return dc

    cpdef double distance_DC(self, double z):
        return self.c/self.H0 * self.__dc(z)

    cpdef double __dz(self, double z):
        cdef double fK
        fK = sqrt(fabs(self.omk))
        if(self.omk > 1e-8):
            return 1/fK * sinh(fK * self.__dc(z))
        elif(self.omk < -1e-8):
            return 1/fK * sin(fK * self.__dc(z))
        
        return self.__dc(z)

    cpdef double distance_DM(self, double z):
        return self.c/self.H0 * self.__dz(z)

    cpdef double distance_DA(self, double z):
        return 1/(1+z)*self.c/self.H0 * self.__dz(z)

    cpdef double distance_DL(self, double z):
        return (1+z) * self.c/self.H0 * self.__dz(z)

    cpdef double distance_DV(self, double z):
        cdef double dv3
        dv3 = pow(self.distance_DM(z), 2) * self.c*z/self.Hofz(z)
        return pow(dv3, 1.0/3.0)

