#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: gapp_kernel.pyx
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2018-12-25 21:29:48
#==================================

cimport cython

import numpy as np
import scipy.optimize as opt
from scipy import integrate
from gapp import dgp
import threading

#@cython.boundscheck(False)
#@cython.wraparound(False)

cdef extern from "math.h":
    double sqrt(double xx)
    double exp(double xx)
    double sinh(double xx)
    double sin(double xx)
    double cos(double xx)
    double fabs(double xx)
    double log10(double xx)
    double log(double xx)

#****************************************************
# give the model parameter use gapp
#****************************************************
cdef class recfuncpar:
    cdef double sigf, lenf

    def __cinit__(self, filename, double[:] initheta):
        '''
        rec_func_par:
        ==========================
        input: filename of data, initial theta
        return: sigmaf and length
        '''
        cdef double[:,:] dats
        
        dats = np.loadtxt(filename, unpack=True)
        self.sigf, self.lenf = self.__output_gp_pars(dats, initheta)

    cpdef double[:] __output_gp_pars(self,double[:,:] dats, double[:] initheta):
        '''
        cpdef double[:] output_gp_pars(double[:,:] dats):
        =================================================
        return theta
        '''
        cdef double zi, zx
        cdef int ns
        cdef double[:] theta
        cdef double[:,:] recH
        
        theta = np.zeros(2,dtype=np.float64)
        
        ns = 2
        zi = dats[0,0]
        zx = dats[0,-1]
        recH = np.zeros((ns,3),dtype=np.float64)
        
        gpr = dgp.DGaussianProcess(dats[0],dats[1],dats[2],cXstar=(zi,zx,ns),\
                theta=initheta)
        recH, theta = gpr.gp()
        return theta

    cpdef double[:] theta(self):
        '''
        theta = [sigf, lenf]
        '''
        return np.array([self.sigf, self.lenf])

#*****************************************************
# gaussian process modules
#*****************************************************
cdef class gaussian_process:
    cdef double[:,:] dats
    cdef double sigmaf, length
    cdef double[:,:] covdat, cov, invcov
    cdef int num    
    cdef double c_light

    def __cinit__(self, filename, covfile=None):
        '''
        self defined Gaussian Process
        '''
        cdef int i, j, cont

        self.dats   = np.loadtxt(filename, unpack=True)
        self.num    = np.shape(self.dats)[1]
        self.covdat = np.zeros((self.num,self.num))
        if(covfile is not None):
            datcov = np.loadtxt(covfile, unpack=True)
            n = datcov.shape[0]
            if(n == self.num):
                for i from 0 <= i < n:
                    self.covdat[i,:] = datcov[i,:]
            elif(n == self.num*self.num):
                cont = 0
                for i from 0<= i<n:
                    for j from 0<=j<n:
                        self.covdat[i,j] = datcov[cont]
                        cont += 1
            else:
                print('too many cov data ...')

        self.c_light = 2.99792458e5 # [km/s]

    def setup_gp(self, double sigf, double lenf):
        self.sigmaf = sigf
        self.length = lenf
        self.cov    = self.__Mcov()
        self.invcov = self.InvMat()

    #***************************
    # about kernel
    #
    cpdef double __kernel(self, double x, double y):
        '''
        Gaussian kernel: sig**2 exp(-(x-y)**2/2/l**2)
        '''
        return self.sigmaf**2*exp(-(x-y)**2/2/self.length**2)
    
    cpdef double __dkerneldx(self, double x, double y):
        '''
        dk(x,y)/dx
        '''
        return (y-x)/self.length**2*self.__kernel(x,y)
    
    cpdef double __dkerneldy(self,double x, double y):
        '''
        dk(x,y)/dy
        '''
        return (x-y)/self.length**2*self.__kernel(x,y)
    
    cpdef double __d2kerneldxdy(self,double x, double y):
        '''
        d^2k(x,y)/dxdy
        '''
        return (1-(x-y)**2/self.length**2)/self.length**2*self.__kernel(x,y)
   
    #******************************
    # about matrix
    #
    cpdef double[:,:] __Mcov(self):
        '''
        M: K(X,X)+C
        Mov(dats,sigf,lenf)
        dats: the datas
        '''
        cdef double[:,:] cov
        cdef double diag
        cdef int i, j, n
        n = self.num
        cov = np.zeros((n,n))
        for i from 0<= i < n:
            for j from 0<= j < n:
                cov[i,j] = self.__kernel(self.dats[0,i],self.dats[0,j]) + self.covdat[i,j]
            diag = cov[i,i] + self.dats[2,i]**2.0
            cov[i,i] = diag
        return cov
    
    cpdef double[:,:] InvMat(self):
        '''
        Inv(M): M = K(X,X) +C
        return M & M^{-1}
        '''
        return np.linalg.inv(self.cov)
    
    #*************************
    # reconstruct func
    #
    cpdef double meanf(self,double x):
        '''
        mean: x = x^*, dats=dats, invcov = InvMat
        return mean
        '''
        cdef int num, i
        cdef double[:] vec
        
        num = self.num
        vec = np.zeros(num)
        for i from 0<= i <num:
            vec[i] = self.__kernel(x, self.dats[0,i])
        return np.dot( np.matmul(vec, self.invcov), self.dats[1,:])
    
    cpdef double dmeanf(self,double x):
        cdef int i, num
        cdef double[:] vec
        
        num = self.num
        vec = np.zeros(num)
        for i in range(num):
            vec[i] = self.__dkerneldx(x, self.dats[0,i] )
        return np.dot( np.matmul(vec, self.invcov), self.dats[1,:])
    
    cpdef double covariance(self,double x, double y):
        '''
        covariance(x,y): cov(fx,fy)
        '''
        cdef int i,num
        cdef double[:] vec1, vec2
        cdef double sig2
        num = self.num
        vec1 = np.zeros(num)
        vec2 = np.zeros(num)
        for i from 0<= i <num:
            vec1[i] = self.__kernel(x,self.dats[0,i])
            vec2[i] = self.__kernel(self.dats[0,i],y)
        sig2 = self.__kernel(x,y) - np.dot(np.matmul(vec1,self.invcov),vec2)
        return sig2

    cpdef double covariance_dfdf(self,double x, double y):
        '''
        covariance_dfdf: cov(df,df)
        '''
        cdef int num, i
        cdef double[:] vec1, vec2
        cdef double sig2
        
        num = self.num
        vec1 = np.zeros(num)
        vec2 = np.zeros(num)
        for i from 0<= i <num:
            vec1[i] = self.__dkerneldx(x, self.dats[0,i])
            vec2[i] = self.__dkerneldy(self.dats[0,i], y)
        sig2 = self.__d2kerneldxdy(x,y) -np.dot(np.matmul(vec1, self.invcov), vec2)
        return sig2
    
    cpdef double covariance_fdf(self,double x, double y):
        '''
        covariancefdf: cov(f,df)
        '''
        cdef int num, i
        cdef double[:] vec1, vec2
        cdef double sig2
        num = self.num
        vec1 = np.zeros(num)
        vec2 = np.zeros(num)
        for i from 0<= i <num:
            vec1[i] = self.__kernel(x, self.dats[0,i])
            vec2[i] = self.__dkerneldx(y, self.dats[0,i])
        sig2 = self.__dkerneldy(x,y) - np.dot( np.matmul(vec1, self.invcov), vec2)
        return sig2

    #****************************
    # return array
    #
    cpdef double[:] array_func(aimfunc, double[:] arr_x):
        '''
        array_func(aimfunc, arr_x)
        ====================================
        return array
        '''
        cdef int n, i
        cdef double[:] arr_func

        n = np.shape(arr_x)[0]
        arr_func = np.zeros(n, dtype=np.float64)

        for i from 0 <= i < n:
            arr_func[i] = aimfunc(arr_x[i])
        return arr_func


    #**************************
    # likelihood
    #
    cpdef double lnlikelihood(self, double[:] theta, out=False):
        '''
        lnlikelihood(double[:] theta):
        '''
        cdef double detcov, lnlike
        cdef double sigf, lenf,chi2, logdet

        sigf, lenf = theta

        self.setup_gp(sigf, lenf)

        detcov = np.linalg.det(self.cov)
        logdet = log(np.abs(detcov))
        chi2 = np.dot(np.matmul(self.dats[1,:],self.invcov),self.dats[1,:])
        lnlike = -0.5* chi2 - 0.5*logdet -self.num/2*log(2*np.pi)
        if(out):
            print('%18.8f'%chi2, '%18.8f'%logdet)
        return lnlike
    
    cpdef double lnprior(self, double[:] theta, double[:] sigma):
        '''
        lnprior
        '''
        cdef double sigf, lenf
        
        sigf, lenf = theta

        if(0<sigf<sigma[0] and 0<lenf<sigma[1]):
            return 0.0
        return -1.0e30
    
    #***************************************
    # lnprobability
    #
    cpdef double lnprobability(self, double[:] theta, double[:] sigma):
        '''
        lnprobability(double[:] theta, double[:] sigma):
        '''
        cdef double lp
        lp = self.lnprior(theta,sigma)
        if(not np.isfinite(lp)):
            return -np.inf
        return self.lnlikelihood(theta)+lp
        
    #******************************************
    # reconstract cosmology background
    #
    cpdef double[:,:] rec_Hz(self, double[:] z):
        '''
        rec_Hz(z):
        =================
        return arr_Hz
        '''
        cdef double[:,:] arr_Hz
        cdef int n, i
        
        n = np.shape(z)[0]
        arr_Hz = np.zeros((3,n),dtype=np.float64)
        
        for i from 0<= i < n:
            arr_Hz[0,i] = z[i]
            arr_Hz[1,i] = self.meanf(z[i])
            arr_Hz[2,i] = np.sqrt(self.covariance(z[i],z[i]))
        return arr_Hz
    
    cpdef double[:,:] cov_Hz(self,double[:] arr_z):
        '''
        cov_Hz(arr_z)
        ======================
        return covhz
        '''
        cdef int i,j,n
        cdef double[:,:] covhz

        n = np.shape(arr_z)[0]
        covhz = np.zeros((n,n), dtype=np.float64)
        for i from 0 <= i < n:
            for j from 0 <= j < n:
                covhz[i,j] = self.covariance(arr_z[i],arr_z[j])
            
        return covhz

    cpdef double __invHz(self,double z):
        '''
        invHz(z)
        =====================
        return 1.0/Hz(z)
        '''
        cdef double dtauda
        dtauda = self.c_light/self.meanf(z)
        return dtauda

    cpdef double[:,:] rec_dc(self,double[:] z, double[:,:] recDc):
        '''
        rec_dc(z)
        ================
        return recDc
        '''
        cdef int n, i
        #cdef double[:,:] recDc
        cdef double fz, err

        n = np.shape(z)[0]
        #recDc = np.zeros((3,n),dtype=np.float64)

        for i from 0 <= i< n:
            if(z[i] == 0.0):
                recDc[0,i] = 0.
                recDc[1,i] = 0.
                recDc[2,i] = 0.
            else:
                recDc[0,i] = z[i]
                recDc[1,i], err = integrate.quad(self.__invHz, 0.0, z[i])
                recDc[2,i] = np.sqrt(self.cov_dcij(z[i],z[i]))

        return recDc


    cpdef double __int_HxHy(self,double x, double y):
        '''
        int_HxHy(x,y)
        =====================
        return 1.0/hx**2/hy**2*covhxhy
        '''
        cdef double covhxhy, hx, hy
        cdef double inthxhy

        covhxhy = self.covariance(x,y)
        hx = self.meanf(x)
        hy = self.meanf(y)
        inthxhy = self.c_light**2/hx**2/hy**2*covhxhy
        return inthxhy

    cpdef double cov_dcij(self,double zi, double zj):
        '''
        cov_dcij(zi,zj)
        ========================
        c^2 \int 1/H(x)^2 1/H(y)^2 cov(H(x),H(y)) dxdy
        '''
        cdef double covij, err

        covij, err = integrate.dblquad(self.__int_HxHy, 0,zi, 0, zj)
        return covij

    cpdef double[:,:] cov_dc(self,double[:] arr_z):
        '''
        cov_dc(arr_z)
        =====================
        return covdc
        '''
        cdef int i,j, n
        cdef double[:,:] covdc
        cdef double fzz

        n = np.shape(arr_z)[0]
        covdc = np.zeros((n,n),dtype=np.float64)

        for i from 0<= i < n:
            print('This is the %s-th row:'%i)
            for j from 0<= j <n:
                if(arr_z[i] ==0.0 and arr_z[j] == 0):
                    covdc[i,j] = 0.0
                else:
                    covdc[i,j] = self.cov_dcij(arr_z[i],arr_z[j])
        print('Done')
        return covdc

    cpdef double[:,:] nstd_cov_dc(self,double[:] arr_x, double[:] arr_y, \
            double[:,:] covdc ):
        '''
        cov_dc(arr_x, arr_y)
        =====================
        return covdc
        '''
        cdef int i,j, n, m

        n = np.shape(arr_x)[0]
        m = np.shape(arr_y)[0]

        for i from 0<= i < n:
            for j from 0<= j <m:
                covdc[i,j] = self.cov_dcij(arr_x[i],arr_y[j])

        return covdc


#========================================================
# GaussianProcess

import numpy as np

cdef class GaussianProcess:
    '''
    GaussianProcess(X,Y,Cov)
    =============================
    X: data
    Y: data
    Cov: covariance matrix
    '''
    cdef double sigf, lenf, logDet
    cdef double[:] X, Y
    cdef double[:,:] Cov, Mcov, IMcov
    cdef int n
    
    def __cinit__(self, double[:] X, double[:] Y, double[:,:] Cov):
        self.X = X
        self.Y = Y
        self.Cov = Cov
        self.n = np.shape(X)[0]
    
    def setup_gp(self, double sigf, double lenf):
        self.sigf = sigf
        self.lenf = lenf
        self.__setdata()
    
    def __kernel(self, double x, double y):
        cdef double pdf
        pdf = self.sigf**2*np.exp(-0.5*(x-y)**2/self.lenf**2)
        return pdf
    
    def __setdata(self):
        cdef double[:,:] L, B, Eye
        cdef int i, j
        self.Mcov = np.zeros((self.n, self.n))
        Eye = np.eye(self.n)
        for i from 0<=i<self.n:
            for j from 0<=j<self.n:
                self.Mcov[i,j] = self.Cov[i,j] +self.__kernel(self.X[i], self.X[j])
        L = np.linalg.cholesky(self.Mcov)
        B = np.linalg.solve(L, Eye)
        self.IMcov = np.linalg.solve(np.transpose(L), B)
        self.logDet = 2*np.sum(np.log(np.diagonal(L)))
        
    cpdef double[:,:] __InvMat(self, double[:,:] mat):
        cdef double[:,:] L, B, Eye, imat
        cdef int n
        n = np.shape(mat)[0]
        Eye = np.ones(n)
        L = np.linalg.cholesky(mat)
        B = np.linalg.solve(L, Eye)
        imat = np.linalg.solve(np.transpose(L), B)
        return imat
    
    cpdef double mean(self, double x):
        '''
        mean(x)
        ==========================
        '''
        cdef double[:] vec
        cdef double mu
        vec = np.array([self.__kernel(x, y) for y in self.X])
        mu = np.dot( np.matmul(vec, self.IMcov), self.Y)
        return mu
    
    cpdef double covariance(self, double x, double y):
        '''
        covariance(x,y)
        ==============================
        '''
        cdef double[:] vec1, vec2
        cdef double sig2
        vec1 = np.array([self.__kernel(x, z) for z in self.X])
        vec2 = np.array([self.__kernel(z, y) for z in self.X])
        sig2 = self.__kernel(x,y) - np.dot(np.matmul(vec1, self.IMcov), vec2)
        return sig2        
        
