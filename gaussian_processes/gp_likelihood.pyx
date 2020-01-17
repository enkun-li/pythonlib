#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: gp_likelihood.pyx
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-03-20 20:35:58
#==================================

cimport cython
import numpy as np
import scipy.optimize as opt
from scipy import integrate
import sys
if(not sys.path[0][-10:-1] == 'pythonlib'):
    sys.path.insert(0,'/home/ekli/myworks/pythonlib/')
from integratefunc import qgaus, qqgaus, qgaus2d_fix

cdef extern from "math.h":
    double sqrt(double xx)
    double exp(double xx)
    double sinh(double xx)
    double sin(double xx)
    double cos(double xx)
    double fabs(double xx)
    double log10(double xx)
    double log(double xx)


#===================================================
# gaussian process
#
cdef class GaussianProcesses:
    '''
    GaussianProcesses(X,Y,Err)
    ================================================
    X: data x direction
    Y: data y direction
    Err: nxn covariance matrix
    '''
    cdef double[:] X, Y
    cdef double[:,:] Cov, Mcov, IMcov
    cdef int n
    cdef double sigf, lenf, logDet
    
    def __cinit__(self, double[:] X, double[:] Y, double[:,:] Cov, 
            sigf=None, lenf=None):
        self.n = np.shape(X)[0]
        self.X = X
        self.Y = Y
        self.Cov = Cov
        if(sigf is not None and lenf is not None):
            self.setup_gp(sigf,lenf)
        return
        
    def setup_gp(self, double sigf, double lenf):
        '''
        setup_gp(sigf, lenf)
        ====================================
        sigf: the strength
        lenf: the length
        '''
        cdef int i, j
        cdef double[:,:] L, B, Eye
        self.sigf = sigf
        self.lenf = lenf
        self.Mcov = np.zeros((self.n, self.n))
        for i from 0<=i<self.n:
            for j from 0<=j<self.n:
                self.Mcov[i,j] = self.Cov[i,j] + self.kernel(self.X[i], self.X[j])
        #self.IMcov = Matrix_Inverse_Chol(self.Mcov)
        #self.logDet = Matrix_logDet(self.Mcov)
        L = np.linalg.cholesky(self.Mcov)
        Eye = np.eye(self.n)
        B = np.linalg.solve(L,Eye)
        self.IMcov = np.linalg.solve(np.transpose(L), B)
        self.logDet = 2*np.sum(np.log(np.diagonal(L)))
        return        
      
    @cython.boundscheck(False)
    cpdef double kernel(self,double x, double y):
        '''
        kernel(x,y)
        ======================================
        return: sigf^2 exp(-(x-y)^2/2/lenf^2)
        '''
        cdef double k
        k = self.sigf**2*exp(-(x-y)**2/2/self.lenf**2)
        return k
    
    @cython.boundscheck(False)
    cpdef double log_likelihood(self):
        '''
        log_likelihood()
        =====================================
        '''
        cdef double chisq
        chisq = np.dot(np.transpose(self.Y), np.matmul(self.IMcov, self.Y))
        chisq += self.logDet + self.n*log(2*np.pi)
        chisq = -0.5*chisq
        return chisq
    
    @cython.boundscheck(False)
    cpdef double __dkerneldx(self, double x, double y):
        '''
        __dkerneldx(x,y)
        ===========================
        dk(x,y)/dx
        '''
        cdef double dk
        dk = (y-x)/self.lenf**2 * self.kernel(x,y)
        return dk

    @cython.boundscheck(False)
    cpdef double __dkerneldy(self, double x, double y):
        '''
        __dkerneldy(x,y)
        ===========================
        dk(x,y)/dy
        '''
        cdef double dk
        dk = (x-y)/self.lenf**2 * self.kernel(x,y)
        return dk

    @cython.boundscheck(False)
    cpdef double __d2kerneldxdy(self, double x, double y):
        '''
        __d2kerneldxdy(x,y)
        ================================
        d2k(x,y)/dxdy
        '''
        cdef double d2k
        d2k = (self.lenf**2 -(x-y)**2)/self.lenf**4 * self.kernel(x,y)
        return d2k

    @cython.boundscheck(False)
    cpdef double __d2kerneldxdx(self, double x, double y):
        '''
        d2kerneldxdx(x,y)
        ================================
        d2k(x,y)/dxdx
        '''
        cdef double d2k
        d2k = -(self.lenf**2 -(x-y)**2)/self.lenf**4 * self.kernel(x,y)
        return d2k
    
    @cython.boundscheck(False)
    cpdef double __d2kerneldydy(self, double x, double y):
        '''
        d2kerneldydy(x,y)
        ================================
        d2k(x,y)/dydy
        '''
        cdef double d2k
        d2k = -(self.lenf**2 -(x-y)**2)/self.lenf**4 * self.kernel(x,y)
        return d2k

    @cython.boundscheck(False)
    cpdef double __d4kerneldx2dy2(self, double x, double y):
        '''
        d4kerneldx2dy2(x,y)
        ================================
        d4k(x,y)/dx2dy2
        '''
        cdef double d4k
        d4k = (3*self.lenf**4 -6*self.lenf**2 *(x-y)**2 +(x-y)**4)/self.lenf**8 \
                * self.kernel(x,y)
        return d4k

    @cython.boundscheck(False)
    cpdef double __diff_kernel(self, double x, double y, int dx=1, int dy=1):
        '''
        __diff_kernel(x,y,dx,dy)
        ==============================
        param x,y: parameter
        param dx,dy: order of derivation
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        n = dx+dy
        d^n k(x,y)/dx^dy dy^dy
        '''
        cdef double diffk
        diffk = 0.0
        if ( dx == 0 and dy == 0):
            diffk = self.kernel(x,y)
        elif(dx == 1 and dy == 0):
            diffk = (y-x)/self.lenf**2 * self.kernel(x,y)
        elif(dx == 2 and dy == 0):
            diffk = -(self.lenf**2 -(x-y)**2)/self.lenf**4 \
                    * self.kernel(x,y)
        elif(dx == 0 and dy == 1):
            diffk = (x-y)/self.lenf**2 * self.kernel(x,y)
        elif(dx == 1 and dy == 1):
            diffk = (self.lenf**2 -(x-y)**2)/self.lenf**4 \
                    * self.kernel(x,y)
        elif(dx == 2 and dy == 1):
            diffk = -(3*self.lenf**2 -(x-y)**2) *(x-y)/self.lenf**6 \
                    * self.kernel(x,y)
        elif(dx == 0 and dy == 2):
            diffk = -(self.lenf**2 -(x-y)**2)/self.lenf**4 \
                    * self.kernel(x,y)
        elif(dx == 1 and dy == 2):
            diffk = (3*self.lenf**2 -(x-y)**2) *(x-y)/self.lenf**6 \
                    * self.kernel(x,y)
        elif(dx == 2 and dy == 2):
            diffk = (3*self.lenf**4 -6*self.lenf**2 *(x-y)**2 +(x-y)**4) \
                    /self.lenf**8 * self.kernel(x,y)
        else:
            print('There is no higher order of diff kernel.')
            sys.exit(0)
        return diffk

#====================================================
# find the best-fit value of sigf and lenf
#
cdef class min_GaussianProcesses(GaussianProcesses):
    '''
    min_GaussianProcesses(X,Y,Err)
    ==============================
    '''

    def __cinit__(self, double[:] X, double[:] Y, double[:,:] Cov):
        super(GaussianProcesses, self).__init__(X,Y,Cov)
        return

    @cython.boundscheck(False)
    cpdef double prob(self, theta):
        cdef double sigf, lenf, chisq
        sigf,lenf = theta
        self.setup_gp(sigf,lenf)
        chisq = -self.log_likelihood()
        return chisq
    
    def min_loglikelihood(self,sigmax=None,lenmax=None):
        bnds = ((0,sigmax), (0,lenmax))
        theta0 = (np.max(self.Y), np.max(self.X))
        res = opt.minimize(self.prob, theta0, bounds=bnds)
        return res
        
#====================================================
# gaussian process reconstract
#
cdef class rec_GaussianProcesses(GaussianProcesses):
    '''
    rec_GaussianProcesses(X,Y,Err)
    ==============================
    '''

    def __cinit__(self, double[:] X, double[:] Y, double[:,:] Cov, double sigf, double lenf):
        super(GaussianProcesses, self).__init__(X,Y,Cov,sigf,lenf)
        #self.setup_gp(sigf, lenf)
        return

    def rec_mu(self,x):
        '''
        rec_mu(x)
        ===========================
        param x: a single value or array
        '''
        cdef int i, j, k, n

        try:
            n = len(x)
        except:
            n = 1
        
        if(n==1):
            mu = 0.0
            for i from 0<= i < self.n:
                for j from 0<= j < self.n:
                    mu += self.kernel(x, self.X[i]) * self.IMcov[i,j] *self.Y[j]
        else:
            mu = np.zeros(n)
            for k from 0 <= k < n:
                for i from 0 <= i < self.n:
                    for j from 0 <= j <self.n:
                        mu[k] += self.kernel(x[k], self.X[i]) * self.IMcov[i,j] *self.Y[j]
        return mu

    def rec_dmu(self,x):
        '''
        rec_dmu(x)
        ===============================
        param x: a single value or array
        '''
        cdef int i, j, k, n

        try:
            n=len(x)
        except:
            n=1

        if(n==1):
            mu = 0.0
            for i from 0 <= i < self.n:
                for j from 0<= j < self.n:
                    mu += self.__dkerneldx(x, self.X[i]) * self.IMcov[i,j] *self.Y[j]
        else:
            mu = np.zeros(n)
            for k from 0 <= k < n:
                for i from 0 <= i < self.n:
                    for j from 0 <= j <self.n:
                        mu[k] += self.__dkerneldx(x[k], self.X[i]) * self.IMcov[i,j] *self.Y[j]
        return mu
    
    def rec_ddmu(self,x):
        '''
        rec_ddmu(x)
        ===============================
        param x: a single value or array
        '''
        cdef int i, j, k, n

        try:
            n=len(x)
        except:
            n=1

        if(n==1):
            mu = 0.0
            for i from 0 <= i < self.n:
                for j from 0<= j < self.n:
                    mu += self.__d2kerneldxdx(x, self.X[i]) * self.IMcov[i,j] *self.Y[j]
        else:
            mu = np.zeros(n)
            for k from 0 <= k < n:
                for i from 0 <= i < self.n:
                    for j from 0 <= j <self.n:
                        mu[k] += self.__d2kerneldxdx(x[k], self.X[i]) * self.IMcov[i,j] *self.Y[j]
        return mu

    def rec_covarianve(self, x):
        '''
        rec_covarianve(x)
        ===========================
        param x: a single value or array

        '''
        cdef int i, j, k, l

        try:
            n = len(x)
        except:
            n = 1
            
        if(n==1):
            cov = self.kernel(x,x)
            for i from 0 <= i < self.n:
                for j from 0 <= j < self.n:
                    cov -= self.kernel(x, self.X[i]) * self.IMcov[i,j] \
                            *self.kernel(self.X[j], x)
        else:
            cov = np.zeros((n,n))
            for k from 0 <= k < n:
                for l from 0 <= l < n:
                    cov[k,l] = self.kernel(x[k], x[l])
                    for i from 0 <= i < self.n:
                        for j from 0 <= j < self.n:
                            cov[k,l] -= self.kernel(x[k], self.X[i]) \
                                    *self.IMcov[i,j] * self.kernel(self.X[j], x[l])
        return cov
        
    def rec_dfdfcovarianve(self, x):
        '''
        rec_dfdfcovarianve(x)
        ===========================
        param x: a single value or array

        ''' 
        cdef int i, j, k, l

        try:
            n = len(x)
        except:
            n = 1
            
        if(n==1):
            cov = self.__d2kerneldxdy(x,x)
            for i from 0 <= i < self.n:
                for j from 0 <= j < self.n:
                    cov -= self.__dkerneldx(x, self.X[i]) * self.IMcov[i,j] \
                            *self.__dkerneldy(self.X[j], x)
        else:
            cov = np.zeros((n,n))
            for k from 0 <= k < n:
                for l from 0 <= l < n:
                    cov[k,l] = self.__d2kerneldxdy(x[k], x[l])
                    for i from 0 <= i < self.n:
                        for j from 0 <= j < self.n:
                            cov[k,l] -= self.__dkerneldx(x[k], self.X[i]) \
                                    *self.IMcov[i,j] * self.__dkerneldy(self.X[j], x[l])
        return cov
     
    def rec_ddfddfcovarianve(self, x):
        '''
        rec_ddfddfcovarianve(x)
        ===========================
        param x: a single value or array

        '''
        cdef int i,j,k,l

        try:
            n = len(x)
        except:
            n = 1
            
        if(n==1):
            cov = self.__d4kerneldx2dy2(x,x)
            for i from 0 <= i < self.n:
                for j from 0 <= j < self.n:
                    cov -= self.__d2kerneldxdx(x, self.X[i]) * self.IMcov[i,j] \
                            *self.__d2kerneldydy(self.X[j], x)
        else:
            cov = np.zeros((n,n))
            for k from 0 <= k < n:
                for l from 0 <= l < n:
                    cov[k,l] = self.__d4kerneldx2dy2(x[k], x[l])
                    for i from 0 <= i < self.n:
                        for j from 0 <= j < self.n:
                            cov[k,l] -= self.__d2kerneldxdx(x[k], self.X[i]) \
                                    *self.IMcov[i,j] * self.__d2kerneldydy(self.X[j], x[l])
        return cov   

    def rec_diff_covarianve(self, x, dx=1,dy=1):
        '''
        rec_diff_covarianve(x,dx,dy)
        ===========================
        param x: a single value or array

        '''
        cdef int i,j,k,l

        try:
            n = len(x)
        except:
            n = 1
            
        if(n==1):
            cov = self.__diff_kernel(x,x,dx,dy)
            for i from 0 <= i < self.n:
                for j from 0 <= j < self.n:
                    cov -= self.__diff_kernel(x, self.X[i],dx,0) * self.IMcov[i,j] \
                            *self.__diff_kernel(self.X[j], x,0,dy)
        else:
            cov = np.zeros((n,n))
            for k from 0 <= k < n:
                for l from 0 <= l < n:
                    cov[k,l] = self.__diff_kernel(x[k], x[l], dx, dy)
                    for i from 0 <= i < self.n:
                        for j from 0 <= j < self.n:
                            cov[k,l] -= self.__diff_kernel(x[k], self.X[i], dx, 0) \
                                    *self.IMcov[i,j] * self.__diff_kernel(self.X[j], x[l], 0, dy)
        return cov

    @cython.boundscheck(False)
    cpdef double rec_mu_x(self, double x):
        '''
        rec_mu_x(x)
        ===========================
        param x: a single value
        '''
        cdef double mu
        cdef int i, j

        mu = 0.0
        for i from 0 <= i < self.n:
            for j from 0 <= j < self.n:
                mu += self.kernel(x, self.X[i]) * self.IMcov[i,j] * self.Y[j]
        return mu
    
    @cython.boundscheck(False)
    cpdef double rec_diff_covarianve_x_y(self, double x, double y, int dx=1, int dy=1):
        '''
        rec_diff_covarianve(x,y,dx,dy)
        ===========================
        param x: a single value or array

        '''
        cdef double cov
        cdef int i, j

        cov = 0.0
        for i from 0 <= i < self.n:
            for j from 0 <= j < self.n:
                cov += self.__diff_kernel(x, self.X[i], dx, 0) * self.IMcov[i,j] \
                        * self.__diff_kernel(self.X[j], y, 0, dy)
        cov = self.__diff_kernel(x,y,dx,dy) -cov
        return cov

#====================================================
# reconstract extern functions
#
cdef class rec_extern_GaussianProcesses(rec_GaussianProcesses):
    '''
    rec_extern_GaussianProcesses(rec_GaussianProcesses)
    =========================================================
    '''
    
    def __cinit__(self, double[:] X, double[:] Y, double[:,:] Cov, double sigf, double lenf):
        super(rec_GaussianProcesses, self).__init__(X,Y,Cov,sigf,lenf)
        return

    def integrate_one_over_mu(self, double xe, double xs=0):
        '''
        integrate_one_over_mu(xs, xe)
        ====================================
        param xe: up limit (end)
        param xs: down limit (start)
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ff = \int_xs^xe 1/mu(x) dx
        '''
        cdef double ff, err, x
        func = lambda x: 1/self.rec_mu_x(x)
        #ff, err = integrate.quad(func, xs, xe)
        ff = qgaus(func, xs, xe)
        return ff

    def integrate_covx_over_mu2(self, double xe1, double xe2, double xs=0):
        '''
        integrate_covx_over_mu(xe1, xe2, xs)
        ===================================
        '''
        cdef double ff, err, x, y
        func = lambda x,y: self.rec_diff_covarianve_x_y(x,y,dx=0,dy=0) / \
                self.rec_mu_x(x)**2 / self.rec_mu_x(y)
        ff = qgaus(func, xs, xe1, args=(xe2))
        return ff

    def integrate_db_covxy_over_mu2(self, double xe1, double xe2, double xs=0):
        '''
        integrate_db_covxy_over_mu2(xe1,xe2,xs)
        ===========================================
        '''
        cdef double ff, x, y
        func = lambda x, y: self.rec_diff_covarianve_x_y(x,y,dx=0,dy=0) / \
                self.rec_mu_x(x)**2 / self.rec_mu_x(y)**2
        ff = qgaus2d_fix(func, xs, xe1, xs, xe2)
        return ff

    def Cov_integrate_one_over_mu(self, double xe1, double xe2, double xs):
        '''
        Cov_integrate_one_over_mu(xe1, xe2, xs)
        =====================================
        param xe1: up limit of mu_1
        param xe2: up limit of mu_2
        param xs:  down limit of all
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        covariance for integrate_one_over_mu
        '''
        cdef double ff, err, x, y
        func = lambda y, x: 1/self.rec_mu_x(x)**2/self.rec_mu_x(y)**2 \
                *self.rec_diff_covarianve_x_y(x,y,dx=0,dy=0)
        #ff, err = integrate.dblquad(func, xs, xe1, xs, xe2)
        ff = qgaus2d_fix(func, xs, xe1, xs, xe2)
        return ff

    def integrate_f_mu_x(self, double xe, double xs=0, 
            double alp=-1, double bet=0.0):
        '''
        integrate_f_mu_x(xe,xs,alp,bet)
        ======================================
        param xe: up limit (end)
        param xs: down limit (start)
        param alp: power of mu
        param bet: power of (1+x)
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        \int_xs^xs mu^alp (1+x)^bet dx
        '''
        cdef double ff, err, x
        func = lambda x: self.rec_mu_x(x)**alp * (1+x)**bet
        #ff, err = integrate.quad(func, xs, xe)
        ff = qgaus(func, xs, xe)
        return ff

    def Cov_integrate_f_mu_x(self, double xe1, double xe2, double xs, 
            double alp=-1, double bet=0.0):
        '''
        Cov_integrate_f_mu_x(xe1, xe2, xs)
        =====================================
        param xe1: up limit of mu_1
        param xe2: up limit of mu_2
        param xs:  down limit of all
        param alp: power of mu
        param bet: power of (1+x)
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        covariance for integrate_f_mu_x
        '''
        cdef double ff, err, x, y
        func = lambda y, x: alp*self.rec_mu_x(x)**(alp-1)*(1+x)**bet \
                * alp*self.rec_mu_x(y)**(alp-1) *(1+x)**bet\
                * self.rec_diff_covarianve_x_y(x,y,dx=0,dy=0)
        #ff, err = integrate.dblquad(func, xs, xe1, xs, xe2)
        ff = qgaus2d_fix(func, xs, xe1, xs, xe2)
        return ff
    
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    @cython.boundscheck(False)
    cpdef double __func(self, double x, func=None):
        '''
        __func(x, func)
        =============================
        '''
        cdef double ff
        if(func is not None):
            ff = func(x)
        else:
            ff = 1.0
        return ff

    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    @cython.boundscheck(False)
    cpdef double int_func(self, double x, func, double alpha=1):
        '''
        int_func(alpha, func)
        ===========================
        mu^\alpha(x) * func(x)
        '''
        cdef double ff
        ff = self.__func(x, func) * self.rec_mu_x(x)**alpha
        return ff

    @cython.boundscheck(False)
    cpdef double __cov_over_mu(self, double x, double y):
        '''
        __cov_over_mu(x,y)
        ==========================
        Cov(mu_i. mu_j)/mu_i/mu_j
        '''
        cdef double ff
        ff = self.rec_diff_covarianve_x_y(x,y,dx=0,dy=0)/self.rec_mu_x(x)/self.rec_mu_x(y)
        return ff

    @cython.boundscheck(False)
    cpdef double int_funcx(self, double x, double y, func, double alpha=1):
        '''
        int_func(x,y,func,alpha)
        ====================================
        func(x) * H^alpha * Cov(Hx, Hy)
        '''
        cdef double ff
        ff = self.int_func(x, func, alpha) * self.__cov_over_mu(x,y)
        return ff

    @cython.boundscheck(False)
    cpdef double int_funcxy(self, double y, double x, func, double alpha=1):
        '''
        __int_func(y,x,func,alpha)
        ====================================
        func1(x) * H^alpha * func2(x) * H^beta * Cov(Hx, Hy)
        '''
        cdef double ff
        
        ff = self.int_func(x, func, alpha) * self.int_func(y, func, alpha) \
                * self.__cov_over_mu(x,y)
        return ff

    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    @cython.boundscheck(False)
    cpdef double integrate_f_pq_x(self, Afunc, Bfunc, double xe, double xs, 
            double alp=-1, double bet=0.0):
        '''
        integrate_f_pq_x(Afunc,Bfunc,xe,xs,alp,bet)
        ==================================================
        Afunc(x)*mu(x)^bet *\int_0^x Bfunc(x') * mu(x')^alp dx'
        '''
        cdef double ff1, ff2, err, ff

        ff1 = self.int_func(xe, Afunc, bet)
        #ff2, err = integrate.quad(self.int_func, xs, xe, args=(Bfunc, alp))
        ff2 = qgaus(self.int_func, xs, xe, args=(Bfunc, alp))
        ff = ff1 *ff2
        return ff

    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    @cython.boundscheck(False)
    cpdef double Cov_integrate_f_pq_x(self, Afunc, Bfunc, double xe1, double xe2,
            double xs, double alp=1, double bet=0.0):
        '''
        Cov_integrate_f_pq_x(Afunc,Bfunc,xe1,xe2,xs,alp,bet)
        ========================================================
        '''
        cdef double ff1, ff2, ff3, ff4, ff, err
        ff1 = self.integrate_f_pq_x(Afunc, Bfunc, xe1, xs, alp, bet) * \
                self.integrate_f_pq_x(Afunc, Bfunc, xe2, xs, alp, bet) * \
                self.__cov_over_mu(xe1,xe2)

        #ff, err = integrate.quad(self.int_funcx, xs, xe2, args=(xe1,Bfunc,alp) )
        ff = qgaus(self.int_funcx, xs, xe2, args=(xe1,Bfunc,alp) )
        ff2 = self.integrate_f_pq_x(Afunc, Bfunc, xe1, xs, alp, bet) * \
                self.int_func(xe2, Afunc, bet) * ff

        #ff, err = integrate.quad(self.int_funcx, xs, xe1, args=(xe2,Bfunc,alp))
        ff = qgaus(self.int_funcx, xs, xe1, args=(xe2,Bfunc,alp))
        ff3 = self.int_func(xe1, Afunc, bet) *ff * \
                self.integrate_f_pq_x(Afunc, Bfunc, xe2, xs, alp, bet)
                
        #ff, err = integrate.dblquad(self.int_funcxy, xs, xe1, xs, xe2, 
        #        args=(Bfunc, alp))
        ff = qgaus2d_fix(self.int_funcxy, xs, xe1, xs, xe2, 
                args=(Bfunc, alp))
        ff4 = self.int_func(xe1, Afunc, bet) *self.int_func(xe2, Afunc, bet) *ff

        ff = bet*bet * ff1 + alp*bet* ff2 +alp*bet * ff3 +alp*alp * ff4
        return ff

    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    @cython.boundscheck(False)
    cpdef double __func_x_mu(self, double x, double alp, double bet):
        '''
        __func_x_mu(x, alp, bet)
        =====================================
        '''
        cdef double ff

        ff = (1+x)**alp * self.rec_mu_x(x)**bet
        return ff

    @cython.boundscheck(False)
    cpdef double integrate_func_x_mu(self, double xe, double xs=0, double alp=1.0, double bet=0.0):
        '''
        integrate_func_x_mu(xe,xs,alp,bet)
        ==================================
        \int_0^x (1+x)^alp * mu^bet dx
        '''
        cdef double ff, err
        #ff, err = integrate.quad(self.__func_x_mu, xs, xe, args=(alp, bet))
        ff = qgaus(self.__func_x_mu, xs, xe, args=(alp, bet))
        return ff

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    @cython.boundscheck(False)
    cpdef double __func_x_mu_cov(self, double x, double y, double alp, double bet):
        '''
        __func_x_mu_cov(x,y,alp,bet)
        ======================================
        '''
        cdef double ff

        ff = self.__func_x_mu(x,alp,bet) * self.__cov_over_mu(x,y)
        return ff

    @cython.boundscheck(False)
    cpdef double Cov_integrate_funcx(self, double y, double xe, double xs, 
            double alp=1.0, double bet=1.0):
        '''
        Cov_integrate_funcx(y,xe,xs,alp,bet)
        =============================================
        '''
        cdef double ff, err

        #ff, err = integrate.quad(self.__func_x_mu_cov, xs, xe, 
        #        args=(y,alp,bet))
        ff = qgaus(self.__func_x_mu_cov, xs, xe, 
                args=(y,alp,bet))
        return ff
    
    @cython.boundscheck(False)
    cpdef double __func_yx_mu_cov(self, double y, double x, 
            double alp, double bet):
        '''
        __func_x_mu_cov(x,y,alp,bet)
        ======================================
        '''
        cdef double ff

        ff = self.__func_x_mu(x,alp,bet) * self.__func_x_mu(y,alp,bet) *\
                self.__cov_over_mu(x,y)
        return ff
    
    @cython.boundscheck(False)
    cpdef double Cov_integrate_funcxy(self, double xe1, double xe2, double xs, 
            double alp=1.0, double bet=1.0):
        '''
        Cov_integrate_funcx(xe1,xe2,xs,alp,bet)
        =============================================
        '''
        cdef double ff, err

        #ff, err = integrate.dblquad(self.__func_yx_mu_cov, xs, xe1, xs, xe2,
        #        args=(alp,bet))
        ff = qgaus2d_fix(self.__func_yx_mu_cov, xs, xe1, xs, xe2,
                args=(alp,bet))
        return ff

#====================================================
# reconstract special functions
#
#class rec_special_GaussianProcesses(rec_GaussianProcesses):
#    '''
#    rec_special_GaussianProcesses(X,Y,Cov,sigf,lenf)
#    '''
#    
#    def __init__(self, double[:] X, double[:] Y, double[:,:] Cov, double sigf, double lenf):
#        super(rec_GaussianProcesses, self).__init__(X,Y,Cov,sigf,lenf)
#        return
#

