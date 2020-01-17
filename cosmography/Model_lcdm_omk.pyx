#!/home/ekli/.local/anaconda3/bin/python
#-*- coding: utf-8 -*-  
#==================================
# File Name: Model_lcdm_omk.pyx
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-10-22 10:29:22
#==================================

cimport cython
import numpy as np
from scipy.integrate import quad

cdef extern from "math.h":
    double pow(double x, double y)
    double sqrt(double x)
    double sinh(double x)
    double sin(double x)
    double cos(double x)
    double fabs(double x)

cdef:
    int NPOINT = 5
    double PI = 3.14159265359
    double EPS = 1.0e-11

def gauleg(double x1, double x2):
    cdef:
        int m, i, j
        double z1, z, xm, xl, pp, p3, p2, p1
        double[:] x, w
    
    x = np.zeros([NPOINT+1])
    w = np.zeros([NPOINT+1])
    
    m = int( (NPOINT +1)/2 )
    xm = 0.5*(x2+x1)
    xl = 0.5*(x2-x1)
    
    z1 = 1e30
    for i from 1<=i<=m:
        z = cos(PI*(i-0.25)/(NPOINT+0.5))
        pp = 1.0
        while (fabs(z - z1) > EPS):
            p1 = 1.0
            p2 = 0.0
            for j from 1<=j<=NPOINT:
                p3 = p2
                p2 = p1
                p1 = ((2.0*j - 1.0) *z *p2 -(j-1.0) *p3)/j
            pp = NPOINT*(z*p1 -p2)/(z*z -1.0)
            z1 = z
            z = z1 - p1/pp
        x[i] = xm -xl*z
        x[NPOINT+1-i] = xm+xl*z
        w[i] = 2.0*xl/((1.0 -z*z) *pp *pp)
        w[NPOINT+1-i] = w[i]
    return x, w

cpdef double quad_gau(func, double a, double b, args=()):
    cdef:
        double[:] x, w
        double s
        int i

    x, w = gauleg(a, b)
    s = 0
    for i from 1<= i <=NPOINT:
        s += w[i] * func(x[i], *args)
    return s

cdef class Model_lcdm_omk:
    cdef double H0, omegam, omegak
    cdef double omegab, omegac, omegar, omegag, omeganu
    cdef double z_drag, Tcmb
    cdef double omegax
    cdef double c

    def __cinit__(self, double H0, double omm, double omk=0, 
            double omb=0.0462, double Tcmb = 2.7255, 
            double N_eff = 3.046, feedback=1):
        cdef double b1, b2, ommh2, ombh2
        self.H0 = H0
        self.omegab = omb
        self.omegam = omm
        self.omegac = omm-omb
        self.omegak = omk
        self.omegag = 2.38e-5*(Tcmb/2.7)**4/(H0/100.0)**2
        self.omeganu = 1.62e-5*(Tcmb/2.7)**4/(H0/100.0)**2 *(N_eff/3.0)
        self.omegar = self.omegag+self.omeganu
        self.omegax = 1-omm-omk -self.omegar
        self.c = 2.99792458e5 #[km/s]
        ##----------------------------
        ommh2 = omm*(H0/100)**2
        ombh2 = omb*(H0/100)**2
        b1 = 0.313*(ommh2)**(-0.419) * (1+0.607 *(ommh2)**0.674)
        b2 = 0.238*(ommh2)**0.223
        self.z_drag = 1291*ommh2**0.251/(1+0.659*ommh2**0.828) *(1+b1*ombh2**b2)
        ##------------------------------
        self.Tcmb = Tcmb
        if(feedback ==0):
            print("="*70)
            print("z_drag             = %8.3f"%self.z_drag)
            print("D_M(z_drag)        = %8.3f"%self.distance_DM(self.z_drag))
            print("D_A(z_drag)        = %8.3f"%self.distance_DA(self.z_drag))
            print("dsound_da(z_drag)  = %8.3f"%self.dsound_da(1/(1+self.z_drag)))
            print("rs(z_drag)         = %8.3f"%self.fit_rs_drag())
        return

    cpdef double Hofz(self, double z):
        cdef double hz2
        hz2 = self.omegam*(1+z)**3 +self.omegak*(1+z)**2 +self.omegax
        hz2 += self.omegar*(1+z)**4
        return self.H0 * sqrt(hz2)

    cpdef double __dtauda(self, double z):
        '''
        __dtauda = H0/H(z)
        '''
        return self.H0/self.Hofz(z)

    cpdef double __dc(self, double z):
        cdef double dc, err
        #dc, err = quad(self.__dtauda, 0.0, z)
        dc = quad_gau(self.__dtauda, 0.0, z)
        return dc

    cpdef double distance_DC(self, double z):
        return self.c/self.H0 * self.__dc(z)

    cpdef double __dz(self, double z):
        cdef double fK
        fK = sqrt(fabs(self.omegak))
        if(self.omegak > 1e-8):
            return 1/fK * sinh(fK * self.__dc(z))
        elif(self.omegak < -1e-8):
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
        dv3 = (self.distance_DM(z)**2) * self.c*z/self.Hofz(z)
        return dv3**(1.0/3.0)

    cpdef double dsound_da(self, double a):
        cdef double R, cs, z
        R = 31500*(self.Tcmb/2.7)**(-4)*a*self.omegab*(self.H0/100)**2
        cs = self.c/sqrt(3*(1+R))
        z = 1.0/a -1.0
        return 1.0/self.Hofz(z)/a**2*cs

    cpdef double fit_rs_drag(self):
        cdef double a_drag, rs
        a_drag = 1/(1+self.z_drag)
        rs = quad_gau(self.dsound_da, 1e-8, a_drag)
        return rs

cdef class Model_aalcdm_omk:
    cdef double H0, omegam, omegak
    cdef double omegax
    cdef double c

    def __cinit__(self, double H0, double q0, double j0=1):
        self.H0 = H0
        self.omegam = 2.0/3*(q0+j0)
        self.omegak = 1-j0
        self.omegax = 1-self.omegam -self.omegak
        self.c = 2.99792458e5 #[km/s]
        return

    cpdef double Hofz(self, double z):
        cdef double hz2
        hz2 = self.omegam*(1+z)**3 +self.omegak*(1+z)**2 +self.omegax
        return self.H0 * sqrt(hz2)

    cpdef double dtauda(self, double z):
        return self.H0/self.Hofz(z)

    cpdef double __dc(self, double z):
        cdef double dc, err
        #dc, err = quad(self.dtauda, 0.0, z)
        dc = quad_gau(self.dtauda, 0.0, z)
        return dc

    cpdef double distance_DC(self, double z):
        return self.c/self.H0 * self.__dc(z)

    cpdef double __dz(self, double z):
        cdef double fK
        fK = sqrt(fabs(self.omegak))
        if(self.omegak > 1e-8):
            return 1/fK * sinh(fK * self.__dc(z))
        elif(self.omegak < -1e-8):
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
        dv3 = (self.distance_DM(z)**2) * self.c*z/self.Hofz(z)
        return dv3**(1.0/3.0)



