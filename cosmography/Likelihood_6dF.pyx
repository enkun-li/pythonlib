#!/home/ekli/.local/anaconda3/bin/python
#-*- coding: utf-8 -*-  
#==================================
# File Name: Likelihood_6dF.pyx
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-10-22 13:04:29
#==================================

cimport cython
import numpy as np

cdef extern from "math.h":
    double pow(double x, double y)
    

cdef class Likelihood_6dF:
    cdef double zeff, dvrs, err, rs_rescale

    def __cinit__(self, datafile):
        data = np.loadtxt(datafile, unpack=True)
        self.zeff = data[0]
        self.dvrs = data[1]
        self.err = data[2]
        self.rs_rescale = data[3]
        return
    
    cpdef double loglikelihood(self, DV, double rd=147.78):
        cdef double chisq, dv

        dv = DV(self.zeff)

        if(np.isnan(dv) or dv <0):
            return -3e30

        chisq = ((rd*self.rs_rescale /dv - self.dvrs)/self.err)**2.0
        
        #if(np.isnan(chisq)):
        #    return -3e30
        return -0.5*chisq

