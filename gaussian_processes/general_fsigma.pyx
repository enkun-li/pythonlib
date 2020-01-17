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
    double pow(double xx, double yy)

#===================================================
# gaussian process
#
@cython.final
cdef class GaussianProcesses:
    '''
    GaussianProcesses(X,Y,Err)
    ================================================
    X: data x direction
    Y: data y direction
    Err: nxn covariance matrix
    '''
    cdef double[:] X, Y, xx
    cdef double[:,:] Cov, Mcov, IMcov
    cdef int n, nx
    cdef double sigf, lenf, logDet
    cdef double[:] integrate_xmu, ini_mu
    cdef double[:,:] integrate_covxmu, integrate_db_covxmu, ini_cov_mu
    cdef double[:] fs8_theory
    cdef double[:,:] fs8cov_theory
    
    def __cinit__(self, double[:] X, double[:] Y, double[:,:] Cov, double[:] xx,
            sigf=None, lenf=None):
        self.n = np.shape(X)[0]
        self.X = X
        self.Y = Y
        self.Cov = Cov
        self.xx = xx
        self.nx = xx.shape[0]
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
        L = np.linalg.cholesky(self.Mcov)
        Eye = np.eye(self.n)
        B = np.linalg.solve(L,Eye)
        self.IMcov = np.linalg.solve(np.transpose(L), B)
        self.logDet = 2*np.sum(np.log(np.diagonal(L)))
        return        
      
    @cython.boundscheck(False)
    cpdef double kernel(self,double x, double y): #nogil:
        '''
        kernel(x,y)
        ======================================
        return: sigf^2 exp(-(x-y)^2/2/lenf^2)
        '''
        cdef double k
        
        k = pow(self.sigf,2) * exp(-pow((x-y), 2)/2/ pow(self.lenf, 2) )
        #k = sigf**2 * exp( -(x-y)**2/2/lenf**2 )
        return k
    
    @cython.boundscheck(False)
    cpdef double rec_mu_x(self, double x): #nogil:
        '''
        rec_mu_x(x)
        ===========================
        param x: a single value
        '''
        cdef double mu
        cdef int i, j;

        mu = 0.0
        for i from 0<= i < self.n:
            for j from 0 <= j < self.n:
                mu += self.kernel(x, self.X[i]) * self.IMcov[i,j] * self.Y[j]
        return mu

    @cython.boundscheck(False)
    cpdef double rec_mu_over_mu0_x(self, double x): #nogil:
        '''
        rec_mu_over_mu0_x(x)
        ===========================
        param x: a single value
        '''
        cdef double mu, mu0, ff
        cdef int i, j;

        mu = 0.0
        mu0 = 0.0
        for i from 0<= i < self.n:
            for j from 0 <= j < self.n:
                mu += self.kernel(x, self.X[i]) * self.IMcov[i,j] * self.Y[j]
                mu0 += self.kernel(0.0, self.X[i]) * self.IMcov[i,j] * self.Y[j]
        ff = mu/mu0
        return ff

    @cython.boundscheck(False)
    cpdef double rec_covariance_xy(self, double x, double y): #nogil:
        '''
        rec_covariance_xy(x,y)
        ===========================
        param x: a single value or array

        '''
        cdef double v, cov
        cdef double[:] vecl, vecr, iMY
        cdef int i, j, k, l

        v = 0.0
        for i from 0 <= i < self.n:
            for j from 0 <= j < self.n:
                v += self.kernel(x, self.X[i]) * self.IMcov[i,j] * self.kernel(self.X[j], y)
        cov = self.kernel(x,y) -v
        return cov

    @cython.boundscheck(False)
    cpdef double fx_over_mu(self, double x, double gam): #nogil:
        '''
        fx_over_mu(x, gam)
        =========================
        '''
        cdef double ff
        
        ff = pow((1+x), (3*gam -1))/pow(self.rec_mu_x(x), (2*gam))
        #ff = (1+x)**(3*gam -1)/self.rec_mu_x(x)**(2*gam)
        return ff

    @cython.boundscheck(False)
    cpdef double[:] integrate_fx_over_mu(self, double gam): #nogil:
        '''
        integrate_fx_over_mu(gam)
        ================================
        return array
        '''
        cdef int i
        cdef double[:] ff

        ff = np.zeros(self.nx)

        for i from 0<= i < self.nx:
            ff[i] = qgaus(self.fx_over_mu, 0.0, self.xx[i], args=(gam))
        return ff

    @cython.boundscheck(False)
    cpdef double cov_over_mumu(self, double x, double y): #nogil:
        '''
        cov_over_mumu(x,y)
        =====================
        '''
        cdef double ff

        ff = self.rec_covariance_xy(x,y)/self.rec_mu_x(x)/self.rec_mu_x(y)
        return ff
    
    @cython.boundscheck(False)
    cpdef double cov_mu_over_mu0_xy(self, double x, double y): #nogil:
        '''
        cov_over_mumu(x,y)
        =====================
        '''
        cdef double ff

        ff = self.rec_mu_over_mu0_x(x) * self.rec_mu_over_mu0_x(y) * \
                (self.cov_over_mumu(x,y) - self.cov_over_mumu(0,x) - self.cov_over_mumu(0,y) + self.cov_over_mumu(0.0,0.0))
        return ff

    @cython.boundscheck(False)
    cpdef double fx_over_mu_covx(self, double x, double y, double gam): #nogil:
        '''
        fx_over_mu_covx(x,y,gam)
        ==========================
        '''
        cdef double ff

        ff = pow((1+x), (3*gam -1))/pow(self.rec_mu_x(x), (2*gam)) * self.cov_over_mumu(x,y)
        return ff

    @cython.boundscheck(False)
    cpdef double[:,:] integrate_fx_over_mu_covx(self, double gam): #nogil:
        '''
        integrate_fx_over_mu_covx(gam)
        ======================================
        return matrix
        '''
        cdef int i, j
        cdef double[:,:] ff

        ff = np.zeros((self.nx,self.nx))

        for i from 0<= i < self.nx:
            for j from 0 <= j < self.nx:
                ff[i,j] = qgaus(self.fx_over_mu_covx, 0.0, self.xx[j], args=(self.xx[i], gam))
        return ff

    @cython.boundscheck(False)
    cpdef dbfx_over_mu_covxy(self, double x, double y, double gam): #nogil:
        '''
        dbfx_over_mu_covxy(x,y,gam)
        ===============================
        '''
        cdef double ff

        ff = pow(((1+x)*(1+y)), (3*gam-1))/pow((self.rec_mu_x(x) *self.rec_mu_x(y)), (2*gam)) *self.cov_over_mumu(x,y)
        return ff

    @cython.boundscheck(False)
    cpdef double[:,:] integrate_dbfx_over_mu_covxy(self, double gam): #nogil:
        '''
        integrate_dbfx_over_mu_covxy
        '''
        cdef int i, j
        cdef double[:,:] ff
        
        ff = np.zeros((self.nx, self.nx))

        for i from 0 <= i < self.nx:
            for j from 0 <= j <= i:
                ff[i,j] = qgaus2d_fix(self.dbfx_over_mu_covxy, 0.0, self.xx[i], 0.0, self.xx[j], args=(gam))
                ff[j,i] = ff[i,j]
        return ff

    @cython.boundscheck(False)
    cpdef double[:,:] cov_over_mumu_mat(self):
        '''
        cov_over_mumu_mat()
        ==========================
        '''
        cdef int i, j
        cdef double[:,:] ff

        ff = np.zeros((self.nx, self.nx))

        for i from 0 <= i < self.nx:
            for j from 0 <= j <= i:
                ff[i,j] = self.cov_over_mumu(self.xx[i], self.xx[j])
                ff[j,i] = ff[i,j]

        return ff

    @cython.boundscheck(False)
    cpdef double[:,:] cov_mu_mat(self):
        '''
        cov_mu_mat()
        ==========================
        '''
        cdef int i, j
        cdef double[:,:] ff

        ff = np.zeros((self.nx, self.nx))

        for i from 0 <= i < self.nx:
            for j from 0 <= j <= i:
                ff[i,j] = self.rec_covariance_xy(self.xx[i], self.xx[j])
                ff[j,i] = ff[i,j]

        return ff

    @cython.boundscheck(False)
    cpdef double[:,:] cov_mu_over_mu0_mat(self):
        '''
        cov_over_mumu_mat()
        ==========================
        '''
        cdef int i, j
        cdef double[:,:] ff
        cdef double[:] Eff

        Eff = np.zeros(self.nx)
        ff = np.zeros((self.nx, self.nx))

        for i from 0 <= i < self.nx:
            for j from 0 <= j <= i:
                #if(j == 0):
                #    Eff[j] = self.rec_mu_over_mu0_x(self.xx[j])
                #ff[i,j] = self.cov_over_mumu(self.xx[i], self.xx[j]) \
                #        - self.cov_over_mumu(self.xx[i], 0.0) \
                #        - self.cov_over_mumu(0.0, self.xx[j]) \
                #        + self.cov_over_mumu(0.0, 0.0)
                #ff[i,j] = ff[i,j] * Eff[i] * Eff[j]
                ff[i,j] = self.cov_mu_over_mu0_xy(self.xx[i], self.xx[j])
                ff[j,i] = ff[i,j]

        return ff

    @cython.boundscheck(False)
    cpdef double[:] rec_mu_arr(self): #nogil:
        '''
        rec_mu_arr()
        ==============================
        '''
        cdef int i
        cdef double[:] ff

        ff = np.zeros(self.nx)

        for i from 0 <= i < self.nx:
            ff[i] = self.rec_mu_x(self.xx[i])

        return ff

    @cython.boundscheck(False)
    cpdef double[:] rec_mu_over_mu0_arr(self): #nogil:
        '''
        rec_mu_arr()
        ==============================
        '''
        cdef int i
        cdef double[:] ff

        ff = np.zeros(self.nx)

        for i from 0 <= i < self.nx:
            ff[i] = self.rec_mu_over_mu0_x(self.xx[i])

        return ff

    @cython.boundscheck(False)
    def rec_mu_over_mu0(self): #nogil:
        '''
        rec_mu_over_mu0()
        ==============================
        '''
        cdef int i, j
        cdef double[:] ff, ffi
        cdef double ff0
        cdef double[:,:] fcov, fcovij
        cdef double[:] fcov0i
        cdef double fcov00

        ff = self.rec_mu_arr()
        ff0 = self.rec_mu_x(0.0)
        ffi = np.zeros(self.nx)

        fcov = self.cov_over_mumu_mat()
        fcov00 = self.cov_over_mumu(0.0, 0.0)
        fcov0i = np.zeros(self.nx)
        fcovij = np.zeros((self.nx, self.nx))
        for i from 0<= i < self.nx:
            fcov0i[i] = self.cov_over_mumu(0.0, self.xx[i])

        for i from 0 <= i < self.nx:
            ffi[i] = ff[i]/ff0

        for i from 0 <= i < self.nx:
            for j from 0 <= j <= i:
                fcovij[i,j] = ffi[i]*ff[j]* (fcov[i,j] - fcov0i[i] \
                        - ff[j] + fcov00)
                fcovij[j,i] = fcovij[i,j]
        return ffi, fcovij

    @cython.boundscheck(False)
    cpdef initial_theory_fs8(self, double sig8, double omm, double h0, double gam): #nogil:
        '''
        initial_theory_fs8(gam)
        =========================
        '''
        cdef double B0
        cdef int i, j

        B0 = pow(omm*pow(h0, 2), gam)

        self.integrate_xmu = self.integrate_fx_over_mu(gam)
        self.integrate_covxmu = self.integrate_fx_over_mu_covx(gam)
        self.integrate_db_covxmu = self.integrate_dbfx_over_mu_covxy(gam)
        self.ini_cov_mu = self.cov_over_mumu_mat()
        self.ini_mu = self.rec_mu_arr()

        self.fs8_theory = np.zeros(self.nx)
        self.fs8cov_theory = np.zeros((self.nx, self.nx))
        
        for i from 0 <= i < self.nx:
            self.fs8_theory[i] = sig8 * B0 * pow((1+self.xx[i]), 3*gam) /pow(self.ini_mu[i], 2*gam) \
                    *exp( - B0 * self.integrate_xmu[i])

        for i from 0 <= i < self.nx:
            for j from 0 <= j < self.nx:
                self.fs8cov_theory[i,j] = self.fs8_theory[i] * self.fs8_theory[j] \
                        *4 *pow(gam, 2) * ( self.ini_cov_mu[i,j] \
                        - B0*self.integrate_covxmu[i,j] \
                        -B0* self.integrate_covxmu[j,i] \
                        + pow(B0, 2) * self.integrate_db_covxmu[i,j] )
        return
    
    cpdef double[:] return_fs8_theory(self):
        return self.fs8_theory

    cpdef double[:,:] return_fs8cov_theory(self):
        return self.fs8cov_theory
    
    @cython.boundscheck(False)
    cpdef initial_theory_fs8_dat(self, double sig8, double omm, double h0, double gam): #nogil:
        '''
        initial_theory_fs8_dat(gam)
        ==============================
        the covariance of theory will be ignored
        in this case
        '''
        cdef double B0
        cdef int i

        B0 = pow(omm*pow(h0, 2), gam)

        self.integrate_xmu = self.integrate_fx_over_mu(gam)
        self.ini_mu = self.rec_mu_arr()

        self.fs8_theory = np.zeros(self.nx)
        
        for i from 0 <= i < self.nx:
            self.fs8_theory[i] = sig8 * B0 * pow((1+self.xx[i]), 3*gam) /pow(self.ini_mu[i], 2*gam) \
                    *exp( - B0 * self.integrate_xmu[i])
        return

    @cython.boundscheck(False)
    cpdef double loglikelihood(self, double sig8, double omm, double h0, double gam, 
            double[:] fs8_obs, double[:,:] fs8cov_obs):
        '''
        loglikelihood(sig8, omm, h0, gam, fs8_obs, fs8cov_obs)
        ===========================================================
        '''
        cdef int i, j
        cdef double[:,:] L, B, Eye
        cdef double[:,:] Mcov
        cdef double chisq, chi

        self.initial_theory_fs8(sig8,omm,h0,gam)

        Mcov = fs8cov_obs

        for i from 0<=i<self.n:
            for j from 0<=j<self.n:
                Mcov[i,j] += self.fs8cov_theory[i,j]

        L = np.linalg.cholesky(Mcov)
        Eye = np.eye(self.nx)
        B = np.linalg.solve(L,Eye)
        Mcov = np.linalg.solve(np.transpose(L), B)
        
        chisq  = 0.0
        for i from 0 <= i < self.nx:
            for j from 0 <= j < self.nx:
                chi = (fs8_obs[i] - self.fs8_theory[i]) * Mcov[i,j] \
                        * (fs8_obs[j] - self.fs8_theory[j])
                if(np.isnan(chi) or chi > 1e30):
                    chi = -1e30
                    return chi
                else:
                    chisq += chi

        chisq = -0.5 * chisq
        return chisq

    @cython.boundscheck(False)
    cpdef double loglikelihood_simp(self, double sig8, double omm, double h0, double gam, 
            double[:] fs8_obs, double[:,:] fs8cov_obs):
        '''
        loglikelihood_simp(sig8, omm, h0, gam, fs8_obs, fs8cov_obs)
        ===========================================================
        '''
        cdef int i, j
        cdef double[:,:] L, B, Eye
        cdef double[:,:] Mcov
        cdef double chisq, chi

        self.initial_theory_fs8_dat(sig8,omm,h0,gam)

        Mcov = fs8cov_obs

        L = np.linalg.cholesky(Mcov)
        Eye = np.eye(self.nx)
        B = np.linalg.solve(L,Eye)
        Mcov = np.linalg.solve(np.transpose(L), B)

        chisq  = 0.0
        for i from 0 <= i < self.nx:
            for j from 0 <= j < self.nx:
                chi = (fs8_obs[i] - self.fs8_theory[i]) * Mcov[i,j] \
                        * (fs8_obs[j] - self.fs8_theory[j])
                if(np.isnan(chi) or chi > 1e30):
                    chi = -1e30
                    return chi
                else:
                    chisq += chi

        chisq = -0.5 * chisq
        return chisq
