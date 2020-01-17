#!/home/ekli/.local/anaconda3/bin/python
#-*- coding: utf-8 -*-  
#==================================
# File Name: Likelihood_OHD.pyx
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-10-22 10:13:39
#==================================

cimport cython
import numpy as np

cdef extern from "math.h":
    double pow(double x, double y)
    

cdef class Likelihood_OHD:
    cdef double[:] zeff, Hobs, err
    cdef int nobs

    def __cinit__(self, datafile):
        data = np.loadtxt(datafile, unpack=True)
        self.nobs = np.shape(data)[1]
        self.zeff = data[0]
        self.Hobs = data[1]
        self.err = data[2]
        return
    
    cpdef double loglikelihood(self, Hofz):
        cdef double chisq, dz
        cdef int i
        chisq = 0
        for i from 0<= i < self.nobs:
            dz = Hofz(self.zeff[i])
            if(np.isnan(dz) or dz <0):
                return -3e30

            chisq += ((dz - self.Hobs[i])/self.err[i])**2.0
        
        if(np.isnan(chisq)):
            return -3e30
        return -0.5*chisq

