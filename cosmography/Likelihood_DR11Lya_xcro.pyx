#!/home/ekli/.local/anaconda3/bin/python
#-*- coding: utf-8 -*-  
#==================================
# File Name: Likelihood_DR11Lya_xcro.pyx
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-10-29 11:04:44
#==================================

cimport cython
import numpy as np
import sys
from scipy.interpolate import RectBivariateSpline

cdef extern from "math.h":
    double pow(double x, double y)
    

class Likelihood_DR11_Lya_xcro:
    '''
    Font-Ribeira et al.
    arXiv: 1311.1767
    ===================================
    Data: not Gaussian
        D_H/r_sd,  D_A/r_sd
        mean: 9.0  10.8
        Cov:   0.09     -0.0468
              -0.0468   0.16
    ===================================
    # This likelihood is from Cuceu, Farr, Lemos & Font-Ribera 2019 
    # (1906.11628) and is based on data from Font-Ribera et al. 2014 
    # (1311.1767).
    # Written by James Farr (June 2019).
    BOSS_DR11_Lya_cross.data_directory      = data.path['data']
    BOSS_DR11_Lya_cross.file                = 'BOSS_DR11_scans/BOSS_DR11_scan_aux.txt'
    BOSS_DR11_Lya_cross.xcf_scan            = 'BOSS_DR11_scans/BOSS_DR11_Lya_cross_scan.dat'
    BOSS_DR11_Lya_cross.rd_rescale          = 1.0    
    BOSS_DR11_Lya_cross.transverse_fid      = 11.59
    BOSS_DR11_Lya_cross.parallel_fid        = 8.708
    '''

    def __init__(self):
        datafile = '/home/pub/ekli_work/cosmodata/BAO/BOSS_DR11_scans/BOSS_DR11_Lya_cross_scan.dat'                
        scan = np.loadtxt(datafile, unpack=True)
        
        #Column numbers in scan for data points.
        ap_index = 0
        at_index = 1
        chi2_index = 2
        
        #Get the alphas and make the scan grid.
        ap = np.array(sorted(set(scan[ap_index])))
        at = np.array(sorted(set(scan[at_index])))
        N_ap = ap.shape[0]
        N_at = at.shape[0]
        
        #Add the chi2 column to the grid.
        #Note that the grid is of shape (N_at,N_ap)
        grid = np.reshape(scan[chi2_index], newshape=(N_at,N_ap))
        
        #Make the interpolator (x refers to at, y refers to ap).
        interp = RectBivariateSpline(at,ap,grid,kx=1,ky=1)

        ##Add the dictionary to the object.
        self.interp = interp
        self.zeff = 2.36
        self.transverse_fid = 11.59
        self.parallel_fid = 8.708

        return


    def loglikelihood(self, Hofz, DA, double rd=147.78):
        #cdef double chisq, trans, para
        #cdef int i
        #cdef double 
        c = 2.99792458e5 #[km/s]

        #chisq = 0

        para = c/(Hofz(self.zeff)*rd) / self.parallel_fid
        if(np.isnan(para) or para<0.8 or para > 1.3):
            return -3e30

        trans = DA(self.zeff)/rd / self.transverse_fid
        if(np.isnan(trans) or trans <0.599999 or trans >1.5 ):
            return -3e30

        chisq = self.interp(trans, para)[0,0]
        
        return -0.5*chisq

cdef class Likelihood_DR11_Lya_xcro_gauss:
    '''
    Font-Ribeira et al.
    arXiv: 1311.1767
    ===================================
    Data: not Gaussian
        D_H/r_sd,  D_A/r_sd
        mean: 9.0  10.8
        Cov:   0.09     -0.0468
              -0.0468   0.16
    ===================================
    # This likelihood is from Cuceu, Farr, Lemos & Font-Ribera 2019 
    # (1906.11628) and is based on data from Font-Ribera et al. 2014 
    # (1311.1767).
    # Written by James Farr (June 2019).
    BOSS_DR11_Lya_cross.data_directory      = data.path['data']
    BOSS_DR11_Lya_cross.file                = 'BOSS_DR11_scans/BOSS_DR11_scan_aux.txt'
    BOSS_DR11_Lya_cross.xcf_scan            = 'BOSS_DR11_scans/BOSS_DR11_Lya_cross_scan.dat'
    BOSS_DR11_Lya_cross.rd_rescale          = 1.0    
    BOSS_DR11_Lya_cross.transverse_fid      = 11.59
    BOSS_DR11_Lya_cross.parallel_fid        = 8.708
    '''
    cdef double zeff
    cdef double[:] obs
    cdef double[:,:] Mcov

    def __cinit__(self):
        self.zeff = 2.36
        self.obs = np.array([9.0, 10.8])
        cov = np.array([[0.09, -0.0468], [-0.0468, 0.16]])
        self.Mcov = np.linalg.inv(cov)

        return


    cpdef double loglikelihood(self, Hofz, DA, double rd=147.78):
        cdef double chisq, dh, da
        cdef double[:] vec
        cdef double c = 2.99792458e5 #[km/s]

        dh = c/(Hofz(self.zeff)*rd)
        da = DA(self.zeff)/rd

        if(np.isnan(dh) or dh<0 or np.isnan(da) or da <0):
            return -3e30

        vec = np.array([dh - self.obs[0], da - self.obs[1]])

        chisq = np.dot(vec, np.matmul(self.Mcov, vec))
        
        return -0.5*chisq

