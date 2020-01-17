#!/home/ekli/.local/anaconda3/bin/python
#-*- coding: utf-8 -*-  
#==================================
# File Name: Likelihood_DR11quasars.pyx
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-10-22 14:01:41
#==================================

cimport cython
import numpy as np
import sys

cdef extern from "math.h":
    double pow(double x, double y)
    

cdef class Likelihood_DR11quasars:
    '''
    Font-Ribeira et al.
    arXiv: 1311.1767
    Data: not Gaussian
        D_H/r_sd,  D_A/r_sd
        mean: 9.0  10.8
        Cov:   0.09     -0.0468
              -0.0468   0.16
    '''
    cdef double[:] zeff, obs, err, types
    cdef int nobs

    def __cinit__(self, datafile):
        data = np.loadtxt(datafile, unpack=True)
        self.nobs = data.shape[1]
        self.zeff = data[0]
        self.obs = data[1]
        self.err = data[2]
        self.types = data[3]
        return
    
    cpdef double loglikelihood(self, Hofz, DA, double rd=147.78):
        cdef double chisq, dz
        cdef int i
        cdef double c = 2.99792458e5 #[km/s]

        chisq = 0

        for i from 0<= i <self.nobs:
            if(self.types[i] == 1):
                dz = Hofz(self.zeff[i])
                if(np.isnan(dz) or dz<0):
                    return -3e30
                chisq += ((c/(dz *rd) - self.obs[i])/self.err[i])**2
            elif(self.types[i] == 2):
                dz = DA(self.zeff[i])
                if(np.isnan(dz) or dz<0):
                    return -3e30
                chisq += ((dz/rd -self.obs[i])/self.err[i])**2.0
            else:
                sys.exit("No such type! [Likelihood_DR11quasars]")
        
        #if(np.isnan(chisq)):
        #    return -3e30
        return -0.5*chisq


