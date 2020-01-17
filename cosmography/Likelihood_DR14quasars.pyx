#!/home/ekli/.local/anaconda3/bin/python
#-*- coding: utf-8 -*-  
#==================================
# File Name: Likelihood_DR14quasars.pyx
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-10-22 14:39:46
#==================================

cimport cython
import numpy as np

cdef extern from "math.h":
    double pow(double x, double y)
    

cdef class Likelihood_DR14quasars:
    cdef double zeff, dvrs, err

    def __cinit__(self, datafile):
        data = np.loadtxt(datafile, unpack=True)
        self.zeff = data[0]
        self.dvrs = data[1]
        self.err = data[2]
        return
    
    cpdef double loglikelihood(self, DV, double rd=147.78):
        cdef double chisq, dz
        
        dz = DV(self.zeff)
        if(np.isnan(dz) or dz <0):
            return -3e30    
        chisq = ((dz/rd - self.dvrs)/self.err)**2.0
        
        #if(np.isnan(chisq)):
        #    return -3e30
        return -0.5*chisq

