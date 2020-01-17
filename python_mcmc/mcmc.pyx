#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: mcmc.pyx
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-01-22 21:22:43
#==================================

import numpy as np

cdef class mcmc:
    cdef double[:] pre_theta
    cdef double[:,:] cov
    cdef double pre_lnlike
    cdef int ndim, nsample
    cdef int mode

    def __cinit__(self):
        self.mode = 0
        print('Using MCMC')
        return
    
    def Setup_MCMC(self, aimfunc, double[:,:] initheta, int nsample=3000):
        '''
        Setup_MCMC(aimfunc,initheta,nsample)
        ===============================================================
        aimfunc:  the aim function
        initheta: [mean, min, max, std.err]
        nsample:  the sample number [default: 3000]
        mode:
        0. the general Metroplis Hastings Sampler [default]
        1. Simulated Annealing
        2. estimate of unbiased variance
        '''
        self.ndim = initheta.shape[0]
        self.cov = np.eye(self.ndim)*initheta[:,3]
        self.pre_theta = np.random.uniform(initheta[:,1], initheta[:,2])
        self.pre_lnlike = aimfunc(self.pre_theta)
        self.nsample = nsample
        return
        
    cpdef double[:,:] SimulatedAnnealing(self, aimfunc, int mode=1, \
            double Tmax=100.0, double Tmin=0.01, double Tratio=0.2):
        '''
        SimulateAnnealing(aimfunc,mode,Tmax,Tmin,Tratio)
        =================================================
        aimfunc:  the aim function
        mode:     1-> exp[(pro-pre)/Ti]
                  2-> 1/{1+exp[-(pro-pre)/Ti]}
        Tmax:     the max Temperature [default: 100]
        Tmin:     the min Temperature [default: 0.01]
        Tratio:   the change rate of Temperature [default: 0.2]
        '''
        cdef int i, nt, cont, ni, nj
        cdef double Ti
        cdef double[:,:] result

        self.mode = mode

        nt = -int(np.log(Tmax/Tmin)/np.log(1.0-Tratio)) +1
        cont = nt*self.nsample
        result = np.zeros((self.ndim+1, cont))

        Ti = Tmax
        cont = 0
        while Ti > Tmin:
            ni = cont*self.nsample
            nj = ni+self.nsample
            result[:,ni:nj] = self.MetroplisHastingsSampler(aimfunc,Ti)
            print('Temperature is %5.2f'%Ti)
            Ti = Ti*(1.0-Tratio)
            cont += 1
        return result

    cpdef double[:,:] EstimateVariance(self, aimfunc, double[:] limit):
        '''
        EstimateVariance(aimfunc, limit)
        ===============================================
        aimfunc:
        limit:
        '''
        cdef double mse, MS, yi
        cdef double[:,:] result
        cdef double[:] ybar, sigi
        cdef int i

        
        ybar = np.zeros(self.ndim)
        sigi = np.zeros(self.ndim)

        mse = np.min(limit)
        MS = 10.0
        while MS > mse:
            result = self.MetroplisHastingsSampler(aimfunc)
            for i from 1 <= i < self.ndim:
                ybar[i-1] = np.sum(result[i])/self.nsample
                yi = np.sum(result[i,-500:])/500.0
                sigi[i-1] = (yi-ybar[i-1])**2
            MS = np.max(sigi)
            print('the mini MSE is: %8.4f'%MS)

        return result

    cpdef double[:,:] MetroplisHastingsSampler(self,aimfunc, double Ti=1.0):
        '''
        MetroplisHastingsSampler(aimfunc,Ti)
        ==================================================
        aimfunc:  the function to determin
        Ti:  the Temperature [default 1]
        '''
        cdef double[:] iresult, pro_theta
        cdef double[:,:] result
        cdef double pro_lnlike, alpha, uni
        cdef int i

        iresult = np.zeros(self.ndim+1)
        result = np.zeros((self.ndim+1, self.nsample))

        for i from 0 <= i < self.nsample:
            pro_theta = np.random.multivariate_normal(self.pre_theta, self.cov)
            pro_lnlike = aimfunc(pro_theta)

            alpha = self.__accpetmode(pro_lnlike, self.pre_lnlike, Ti)
            uni = np.random.uniform()
            if(uni < alpha):
                self.pre_theta = pro_theta
                self.pre_lnlike = pro_lnlike
                iresult = np.append(pro_lnlike, pro_theta)
            iresult=np.append(self.pre_lnlike,self.pre_theta)
            result[:,i]=iresult
        return result

    cpdef double __accpetmode(self,double pro, double pre, double Ti):
        '''
        accpetmode(pro,pre,Ti)
        =======================================
        pro:
        pre:
        Ti:
        mode:
        ===> 0: default exp[(pro - pre)]
        ===> 1: exp[(pro-pre)/Ti]
        ===> 2: 1/{1+exp[-(pro-pre)/Ti]}
        '''
        cdef double accept, alpha
        
        if(self.mode==1):
            accept = np.exp((pro-pre)/Ti)
        elif(self.mode==2):
            accept = 1.0/(1.0+np.exp(-(pro-pre)/Ti))
        else:
            accept = np.exp(pro-pre)
        alpha = np.min([1.0,accept])
        return alpha


