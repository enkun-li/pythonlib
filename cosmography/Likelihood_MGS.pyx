#!/home/ekli/.local/anaconda3/bin/python
#-*- coding: utf-8 -*-  
#==================================
# File Name: Likelihood_MGS.pyx
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-10-22 13:19:47
#==================================

cimport cython
import numpy as np

cdef extern from "math.h":
    double pow(double x, double y)
    

cdef class Likelihood_MGS:
    cdef double zeff, dvrs, err
    cdef double[:] prob

    def __cinit__(self, datafile, dataprob):
        data = np.loadtxt(datafile, unpack=True)
        prob = np.loadtxt(dataprob, unpack=True)
        self.zeff = data[0]
        self.dvrs = data[1]
        self.err = data[2]
        self.prob = prob
        return
    
    cpdef double loglikelihood(self, DV, double rd=147.78):
        cdef double rsfidmgs = 148.69, DVfidmgs = 638.9518
        cdef double alpha_min = 0.8005, alpha_max = 1.1985
        cdef int ii
        cdef double alphamgs, chisq, dz
        
        dz = DV(self.zeff)
        if(np.isnan(dz) or dz <0):
            return -3e30

        alphamgs = DV(self.zeff) /rd / (DVfidmgs/rsfidmgs)
        if(alphamgs > alpha_max or alphamgs < alpha_min):
            return -3e30
        else:
            ii = int((alphamgs -alpha_min)/0.001)
            chisq = 0.5*(self.prob[ii] +self.prob[ii+1])

        #if(np.isnan(chisq)):
        #    return -3e30
        return -0.5*chisq
    
    cpdef double loglikelihood_cal(self, DV, double rd=147.78):
        cdef double chisq
        
        chisq = ((DV(self.zeff)/rd - self.dvrs)/self.err)**2.0
        
        if(np.isnan(chisq)):
            return -3e30
        return -0.5*chisq
