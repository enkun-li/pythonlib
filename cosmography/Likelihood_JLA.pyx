#!/home/ekli/.local/anaconda3/bin/python
#-*- coding: utf-8 -*-  
#==================================
# File Name: Likelihood_JLA.pyx
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-10-22 11:02:00
#==================================

cimport cython
import numpy as np

cdef extern from "math.h":
    double pow(double x, double y)
    double log10(double x)
    double log(double x)
    

cdef class Likelihood_JLA:
    cdef double[:] zcmb, zhel, mb, dmb
    cdef double[:,:] Mcov
    cdef int nobs

    def __cinit__(self, datafile, covfile):
        data = np.loadtxt(datafile, unpack=True)
        cov = np.loadtxt(covfile, unpack=True)
        self.nobs = np.shape(data)[1]
        self.zcmb = data[0]
        self.zhel = data[1]
        self.mb = data[2]
        self.dmb = data[3]
        self.Mcov = np.reshape(cov, (self.nobs, self.nobs)) \
                +np.eye(self.nobs)* data[3]**2
        self.Mcov = np.linalg.inv(self.Mcov)
        return
    
    cpdef double loglikelihood(self, DA):
        cdef double chisq, zc, zh
        cdef double[:] difmu
        cdef int i
        chisq = 0

        difmu = np.array([5*log10((1+zc)*(1+zh)*DA(zc)) -mi
            for zc,zh,mi in zip(self.zcmb, self.zhel, self.mb)])
        
        Acov = np.matmul(difmu, self.Mcov)
        magA = np.dot(Acov, difmu)
        magB = np.sum(Acov)
        magC = np.sum(self.Mcov)

        chisq = magA - magB**2/magC # + log(magC/2/np.pi)
        
        if(np.isnan(chisq)):
            return -3e30
        return -0.5*chisq

