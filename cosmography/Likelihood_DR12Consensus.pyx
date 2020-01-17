#!/home/ekli/.local/anaconda3/bin/python
#-*- coding: utf-8 -*-  
#==================================
# File Name: Likelihood_DR12Consensus.pyx
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-10-22 13:38:45
#==================================

cimport cython
import numpy as np
import sys

cdef extern from "math.h":
    double pow(double x, double y)
    

cdef class Likelihood_DR12Consensus:
    cdef double[:] zeff, obs, err
    cdef double[:,:] Mcov
    cdef int nobs

    def __cinit__(self, datafile, datacov):
        data = np.loadtxt(datafile, unpack=True)
        cov = np.loadtxt(datacov, unpack=True)
        self.nobs = np.shape(data)[1]
        self.zeff = data[0]
        self.obs = data[1]
        self.err = data[2]
        self.Mcov = np.linalg.inv(cov)
        return
    
    cpdef double loglikelihood(self, Hofz, DM, double rd=147.78):
        cdef double rs_rescale = .6766815537e-2
        cdef double rs, chisq, dz
        cdef double[:] difbao
        cdef int i

        rs = rd*rs_rescale
        difbao = np.zeros(self.nobs)

        for i from 0<= i < self.nobs:
            if(self.err[i] == 0):
                dz = DM(self.zeff[i])
                if(np.isnan(dz) or dz <0):
                    return -3e30
                difbao[i] = dz/rs - self.obs[i]
            elif(self.err[i] == 1):
                dz = Hofz(self.zeff[i])
                if(np.isnan(dz) or dz <0):
                    return -3e30
                difbao[i] = dz *rs -self.obs[i]
            else:
                sys.exit('There is no such mode Likelihood_DR12Consensus')

        chisq = np.dot( difbao, np.matmul(self.Mcov, difbao) )
        
        #if(np.isnan(chisq)):
        #    return -3e30
        return -0.5*chisq

