#!/home/ekli/.local/anaconda3/bin/python
#-*- coding: utf-8 -*-  
#==================================
# File Name: Likelihood_HST.pyx
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-10-22 10:00:45
#==================================

cimport cython
import numpy as np

cdef extern from "math.h":
    double pow(double x, double y)
    

cdef class Likelihood_HST:
    cdef double zeff, Hobs, err

    def __cinit__(self, datafile):
        data = np.loadtxt(datafile, unpack=True)
        self.zeff = data[0]
        self.Hobs = data[1]
        self.err = data[2]
        return
    
    cpdef double loglikelihood(self, double Hth):
        cdef double chisq
        chisq = ((Hth - self.Hobs)/self.err)**2.0
        return -0.5*chisq



