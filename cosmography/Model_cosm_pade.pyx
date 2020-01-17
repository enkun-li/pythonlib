#!/home/ekli/.local/anaconda3/bin/python
#-*- coding: utf-8 -*-  
#==================================
# File Name: Model_cosm_pade.pyx
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-10-22 16:07:09
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

cdef class Model_cosm_pade:
    cdef double H0, q0, j0, s0, l0, omk
    cdef double omksq
    cdef double c
    cdef double g1,g2,g3,g4,g5,g6,g7,g8,g9
    cdef double f1,f2,f3,f4,f5,f6
    cdef int case #[i: Pi]

    def __cinit__(self, double H0=70.0, double q0=-0.55, double j0=1.0, \
            double s0=-0.35, double l0=3.115, double omk=0, int case=12):
        cdef omksq, omkcb
        self.H0 = H0
        self.q0 = q0
        self.j0 = j0
        self.s0 = s0
        self.l0 = l0
        self.omk = omk
        omksq = omk*omk
        omkcb = omksq *omk
        self.case = case
        self.c = 2.99792458e5 #[km/s]
        #===================
        self.g1 = 1 + q0;
        self.g2 = 2 - j0 + 4*q0 + 3*q0**2 + omk;
        self.g3 = 6 + 18*q0 + 27*q0**2 + 15*q0**3 \
                - j0*(9 + 10*q0) - s0 + 6*(1 + q0)*omk;
        self.g4 = 1 - 2*j0 + 2*q0 + 3*q0**2 + 2*omk;
        self.g5 = 1 + 3*q0 + 8*q0**2 + 6*q0**3 \
                - j0*(5 + 6*q0) - s0 + 2*(1 + q0)*omk;
        self.g6 = 19 + 40*j0**2 - 6*l0 + 76*q0 + 286*q0**2 \
                + 420*q0**3 + 225*q0**4 - 2*j0*(86 + 205*q0 \
                + 150*q0**2) - 66*s0 - 60*q0*s0 - 20*(-2 + j0 \
                - 4*q0 - 3*q0**2)*omk - 14*omksq;
        self.g7 = 2 + 6*q0 + 13*q0**2 + 9*q0**3 - j0*(7 + 8*q0) \
                - s0  + 4*(1 + q0)*omk;
        self.g8 = 2 - 4*j0**2 + 8*q0 + 23*q0**2 + 30*q0**3 \
                + 9*q0**4 - j0*(11 + 25*q0 + 6*q0**2) - 3*s0 \
                - 3*q0*s0 + (2 + 8*j0 + 4*q0 - 6*q0**2)*omk - 4 *omksq;
        self.g9 = 24 + 10*j0**2 - l0 + 96*q0 + 216*q0**2 \
                + 240*q0**3   + 105*q0**4 - j0*(72 + 160*q0 \
                + 105*q0**2) - 16*s0  - 15*q0*s0 + 5*(7 - 2*j0 + 14*q0 \
                + 9*q0**2)*omk + omksq;
        
        self.f1 = 18 + 20*j0**2 - 2*l0 + 72*q0 + 207*q0**2 \
                + 270*q0**3 + 135*q0**4 - j0*(99 + 225*q0 \
                + 160*q0**2) - 27*s0 - 25*q0*s0 - 20*(-2 + j0 \
                - 4*q0 - 3*q0**2)*omk + 2 *omksq;
        self.f2 = 12 - 3*l0 + 60*q0 - 3*l0*q0 + 216*q0**2 \
                + 408*q0**3 + 330*q0**4 + 90*q0**5 \
                - 5*j0**2*(3 + 4*q0) - 38*s0 - 73*q0*s0 \
                - 30*q0**2*s0 - j0*(96 + 326*q0 + 325*q0**2 \
                + 90*q0**3 + 5*s0) + 5*(3 + 9*q0 - 6*q0**3 \
                + j0*(9 + 10*q0) + s0)*omk - 27*(1 + q0) *omksq;
        self.f3 = 12 - 40*j0**3 - 8*l0 + 72*q0 - 16*l0*q0 \
                + 312*q0**2 - 12*l0*q0**2 + 768*q0**3 \
                + 927*q0**4 + 510*q0**5 + 135*q0**6 \
                + j0**2*(-37 - 100*q0 + 40*q0**2) - 68*s0 \
                - 196*q0*s0 - 162*q0**2*s0 - 30*q0**3*s0 \
                - 5*(s0)**2 + 2*j0*(-66 + 2*l0 - 298*q0 - 449*q0**2 \
                - 255*q0**3 - 90*q0**4 - 13*s0 - 20*q0*s0) \
                + 4*(4 + 20*j0**2 - l0 + 16*q0 + 16*q0**2 \
                + 15*q0**4 + j0*(8 + 15*q0 - 30*q0**2) - s0)*omk \
                - 4*(8 + 11*j0 + 16*q0 - 3*q0**2) *omksq + 4 *omkcb;
        self.f4 = 14 - 6*l0 + 70*q0 - 6*l0*q0 + 277*q0**2 \
                + 551*q0**3 + 465*q0**4 + 135*q0**5 \
                - 10*j0**2*(1 + 2*q0) - 61*s0 - 116*q0*s0 \
                - 45*q0**2*s0 - j0*(137 + 472*q0 + 495*q0**2 \
                + 150*q0**3 + 10*s0) + 10*(2 + 6*q0 + q0**2 \
                - 3*q0**3 + j0*(5 + 6*q0) + s0)*omk - 34*(1 + q0) *omksq;
        self.f5 = 4 - 80*j0**3 - 6*l0 + 24*q0 - 12*l0*q0 \
                + 120*q0**2 - 18*l0*q0**2 + 320*q0**3 \
                + 423*q0**4 + 270*q0**5 + 135*q0**6 \
                + 9*j0**2*(1 + 20*q0**2) - 36*s0 - 102*q0*s0 \
                - 78*q0**2*s0 - 15*(s0)**2 +  6*j0*(-10 + 2*l0 \
                - 78*q0**2 - 55*q0**3 - 45*q0**4 - 3*s0 \
                - 2*q0*(23 + 5*s0)) + 6*(3 + 20*j0**2 - 2*l0 + 12*q0 \
                + 42*q0**2 + 60*q0**3 + 45*q0**4 \
                - 2*j0*(12 + 30*q0 + 35*q0**2) - 12*s0 - 10*q0*s0)*omk \
                + 6*(1 - 2*j0 + 2*q0  + 3*q0**2) *omksq - 28 *omkcb;
        self.f6 = 34 + 40*j0**2 - 6*l0 + 136*q0 + 451*q0**2 \
                + 630*q0**3 + 315*q0**4 - j0*(247 + 575*q0 \
                + 390*q0**2) - 81*s0 - 75*q0*s0 + 10*(7 - 2*j0 \
                + 14*q0 + 9*q0**2)*omk - 14 *omksq;

        return

    ## ===================================    
    ## Pade approx of d(z)
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cpdef double pade_cg_dz11(self, double z):
        return z/(1+0.5*self.g1*z);
    
    cpdef double pade_cg_dz12(self, double z):        
        return z/(1 + 0.5*self.g1*z -1.0/12*self.g4*z**2 );
    
    cpdef double pade_cg_dz13(self, double z):
        return z/(1 + 0.5*self.g1 *z -1.0/12*self.g4*z**2 + 1.0/24*self.g5*z**3)

    cpdef double pade_cg_dz14(self, double z):
        return z/(1 + 0.5*self.g1 *z -1.0/12*self.g4*z**2 + 1.0/24*self.g5*z**3 \
            -1.0/720*self.g6*z**4);

    cpdef double pade_cg_dz21(self, double z):
        return (z+self.g4/(6*self.g1)*z**2) / (1+self.g2/(3*self.g1)*z);
    
    cpdef double pade_cg_dz22(self, double z):
        return (z + (self.g5 * z**2)/(2 * self.g4))/(1 + (self.g7 * z)/(2 * self.g4) \
                + (self.g8 * z**2)/(12 * self.g4));
    
    cpdef double pade_cg_dz23(self, double z):
        return (z + (self.g6 * z**2)/(30 * self.g5)) \
            /(1 + (self.f6 * z)/(30 * self.g5) + (self.f4 * z**2)/(60 * self.g5) \
            - (self.f5 * z**3)/(360 * self.g5));
    
    cpdef double pade_cg_dz31(self, double z):
        return (z + (self.g7 * z**2)/(4 * self.g2) \
                - (self.g8 * z**3)/(24 * self.g2)) \
                /(1 + (self.g3 * z)/(4 * self.g2));
    
    cpdef double pade_cg_dz32(self, double z):
        return (z + (self.f4*z**2)/(10*self.g8) + (self.f5*z**3)/(60*self.g8)) \
                /(1 + (self.f2*z)/(5*self.g8) + (self.f3*z**2)/(20*self.g8));
    
    cpdef double pade_cg_dz41(self, double z):
        return (z + (self.f1*z**2)/(10*self.g3) - (self.f2*z**3)/(30*self.g3) \
                + (self.f3*z**4)/(120*self.g3))/(1 + (self.g9*z)/(5*self.g3));
    
    ## *********************************
    ## Pade approx of \pd d(z)/ \pd z
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cpdef double pade_cg_pd_dz11(self, double z):
        return 1/(1+0.5*self.g1*z) - z*0.5*self.g1/(1+0.5*self.g1*z)**2;
    
    cpdef double pade_cg_pd_dz12(self, double z):
        return 1/(1 + 0.5*self.g1*z -1.0/12*self.g4*z**2 ) \
                -z*(0.5*self.g1 -1.0/6*self.g4*z) \
                /(1 + 0.5*self.g1*z -1.0/12*self.g4*z**2)**2;
    
    cpdef double pade_cg_pd_dz13(self, double z):
        return 1/(1 + 0.5*self.g1 *z -1.0/12*self.g4*z**2 \
                + 1.0/24*self.g5*z**3) -z*(0.5*self.g1 -1.0/6*self.g4*z \
                + 1.0/8*self.g5*z**2) /(1 + 0.5*self.g1 *z \
                -1.0/12*self.g4*z**2 + 1.0/24*self.g5*z**3)**2;
    
    cpdef double pade_cg_pd_dz14(self, double z):
        return 1/(1 + 0.5*self.g1 *z -1.0/12*self.g4*z**2 \
                + 1.0/24*self.g5*z**3 -1.0/720*self.g6*z**4) \
                -z*(0.5*self.g1 -1.0/6*self.g4*z + 1.0/8*self.g5*z**2 \
                -1.0/180*self.g6*z**3) /(1 + 0.5*self.g1 *z \
                -1.0/12*self.g4*z**2 + 1.0/24*self.g5*z**3 \
                -1.0/720*self.g6*z**4)**2;
    
    cpdef double pade_cg_pd_dz21(self, double z):
        return (1+self.g4/(3*self.g1)*z) / (1+self.g2/(3*self.g1)*z) \
            -(z+self.g4/(6*self.g1)*z**2) * (self.g2/(3*self.g1)) \
           /(1+self.g2/(3*self.g1)*z)**2;
    
    cpdef double pade_cg_pd_dz22(self, double z):
        return (1 + (self.g5 * z)/(self.g4))/(1 + (self.g7 * z)/(2 * self.g4) \
                + (self.g8 * z**2)/(12 * self.g4)) \
                -(z + (self.g5 * z**2)/(2 * self.g4)) * ((self.g7)/(2 * self.g4) \
                + (self.g8 * z)/(6 * self.g4)) /(1 + (self.g7 * z)/(2 * self.g4) \
                + (self.g8 * z**2)/(12 * self.g4))**2;
    
    cpdef double pade_cg_pd_dz23(self, double z):
        return (1 + (self.g6 * z)/(15 * self.g5)) /(1 + (self.f6 * z)/(30 * self.g5) \
                + (self.f4 * z**2)/(60 * self.g5) - (self.f5 * z**3)/(360* self.g5)) \
                - (z + (self.g6 * z**2)/(30 * self.g5)) * (self.f6/(30 * self.g5) \
                + (self.f4 * z)/(30 * self.g5) - (self.f5 * z**2)/(120 * self.g5)) \
                /(1 + (self.f6 * z)/(30 * self.g5) + (self.f4 * z**2)/(60 * self.g5) \
                - (self.f5 * z**3)/(360 * self.g5))**2;
    
    cpdef double pade_cg_pd_dz31(self, double z):
        return (1 + (self.g7 * z)/(2 * self.g2) - (self.g8 * z**2)/(8 * self.g2)) \
                /(1 + (self.g3 * z)/(4 * self.g2)) - (z + (self.g7 * z**2)/(4 \
                * self.g2) - (self.g8 * z**3)/(24 * self.g2)) * ((self.g3 )/(4 \
                * self.g2)) /(1 + (self.g3 * z)/(4 * self.g2))**2;
    
    cpdef double pade_cg_pd_dz32(self, double z):
        return (1 + (self.f4*z)/(5*self.g8) + (self.f5*z**2)/(20*self.g8)) \
                /(1 + (self.f2*z)/(5*self.g8) + (self.f3*z**2)/(20*self.g8)) \
                - (z + (self.f4*z**2)/(10*self.g8) + (self.f5*z**3)/(60*self.g8)) \
                *((self.f2)/(5*self.g8) + (self.f3*z)/(10*self.g8)) /(1 + (self.f2*z)/(5*self.g8) \
                + (self.f3*z**2)/(20*self.g8))**2;
    
    cpdef double pade_cg_pd_dz41(self, double z):
        return (1 + (self.f1*z)/(5*self.g3) - (self.f2*z**2)/(10*self.g3) \
                + (self.f3*z**3)/(30*self.g3))/(1 + (self.g9*z)/(5*self.g3)) \
                -(z + (self.f1*z**2)/(10*self.g3) - (self.f2*z**3)/(30*self.g3) \
                + (self.f3*z**4)/(120*self.g3)) * ((self.g9)/(5*self.g3)) \
                /(1 + (self.g9*z)/(5*self.g3))**2;
   
    ## ===================================
    cpdef double Hofz(self, double z):        
        cdef double hz
        hz = 0
        if(self.case == 11):
            hz = self.H0/self.pade_cg_pd_dz11(z) \
                    * sqrt(1+self.omk*(self.pade_cg_dz11(z))**2);
        elif(self.case == 12):
            hz = self.H0/self.pade_cg_pd_dz12(z) \
                    * sqrt(1+self.omk*(self.pade_cg_dz12(z))**2);
        elif(self.case == 13):
            hz = self.H0/self.pade_cg_pd_dz13(z) \
                    * sqrt(1+self.omk*(self.pade_cg_dz13(z))**2);
        elif(self.case == 14):
            hz = self.H0/self.pade_cg_pd_dz14(z) \
                    * sqrt(1+self.omk*(self.pade_cg_dz14(z))**2);
        elif(self.case == 21):
            hz = self.H0/self.pade_cg_pd_dz21(z) \
                    * sqrt(1+self.omk*(self.pade_cg_dz21(z))**2);
        elif(self.case == 22):
            hz = self.H0/self.pade_cg_pd_dz22(z) \
                    * sqrt(1+self.omk*(self.pade_cg_dz22(z))**2);
        elif(self.case == 23):
            hz = self.H0/self.pade_cg_pd_dz23(z) \
                    * sqrt(1+self.omk*(self.pade_cg_dz23(z))**2);
        elif(self.case == 31):
            hz = self.H0/self.pade_cg_pd_dz31(z) \
                    * sqrt(1+self.omk*(self.pade_cg_dz31(z))**2);
        elif(self.case == 32):
            hz = self.H0/self.pade_cg_pd_dz32(z) \
                    * sqrt(1+self.omk*(self.pade_cg_dz32(z))**2);
        elif(self.case == 41):
            hz = self.H0/self.pade_cg_pd_dz41(z) \
                    * sqrt(1+self.omk*(self.pade_cg_dz41(z))**2);
        else:
            sys.exit("No such model %s! [Model_cosm_tay]"%self.case)
        return hz

    cpdef double __dz(self, double z):
        cdef double dz
        dz = 0
        if(self.case == 11):
            dz = self.pade_cg_dz11(z)
        elif(self.case == 12):
            dz = self.pade_cg_dz12(z)
        elif(self.case == 13):
            dz = self.pade_cg_dz13(z)
        elif(self.case == 14):
            dz = self.pade_cg_dz14(z)
        elif(self.case == 21):
            dz = self.pade_cg_dz21(z)
        elif(self.case == 22):
            dz = self.pade_cg_dz22(z)
        elif(self.case == 23):
            dz = self.pade_cg_dz23(z)
        elif(self.case == 31):
            dz = self.pade_cg_dz31(z)
        elif(self.case == 32):
            dz = self.pade_cg_dz32(z)
        elif(self.case == 41):
            dz = self.pade_cg_dz41(z)
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

