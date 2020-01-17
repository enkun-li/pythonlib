#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: Matrix_utils.pyx
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-03-23 08:36:21
#==================================

import numpy as np

cdef extern from "math.h":
    double sqrt(double xx)
    double exp(double xx)
    double sinh(double xx)
    double cosh(double xx)
    double sin(double xx)
    double cos(double xx)
    double fabs(double xx)
    double log10(double xx)
    double log(double xx)

cpdef double[:,:] Matrix_Choleskey(M):
    '''
    Matrix_Choleskey(M)
    =========================
    M: positive definite matrix
    '''
    cdef double[:,:] L
    try:
        L = np.linalg.cholesky(M)
    except np.linalg.LinAlgError:
        L = None
        print('Matrix_Choleskey: M is not positive definite.')
    return L

cpdef double[:,:] Matrix_Inverse_Chol(M):
    '''
    Matrix_Inverse_Chol(M)
    ==============================
    M: positive definite matrix
    ---> M = LL*
    ---> M^-1 = L*^-1 L^-1
    ---> LB = I
    ---> B = L^-1L B = L^-1 I
    ---> B = L* L*^-1L^-1 I = L* alpha
    ---> alpha = M^-1
    '''
    cdef double[:,:] Eye, L, B, invM
    cdef int ndim
    ndim = np.shape(M)[0]
    Eye = np.eye(ndim)
    L = Matrix_Choleskey(M)
    B = np.linalg.solve(L, Eye)
    invM = np.linalg.solve(np.transpose(L), B)
    return invM

cpdef double Matrix_logDet(M):
    '''
    Matrix_logDet(M)
    ===============================
    M: positive definite matrix
    ---> M = L*L
    ---> log(|M|) = sum( log( np.diagonal(L) ) )
    '''
    cdef double[:,:] L
    cdef double[:] T
    cdef double logdet
    L = Matrix_Choleskey(M)
    T = np.diagonal(L)
    logdet = 2*np.sum(np.log(T))
    return logdet




