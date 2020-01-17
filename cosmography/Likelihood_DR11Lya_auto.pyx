#!/home/ekli/.local/anaconda3/bin/python
#-*- coding: utf-8 -*-  
#==================================
# File Name: Likelihood_DR11Lya_auto.pyx
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
    

class Likelihood_DR11_Lya_auto:
    '''
    Font-Ribeira et al.
    arXiv: 1311.1767
    ===================================
    Data: not Gaussian
        D_H/r_sd,  D_A/r_sd
        mean: 9.18, 11.27994012
        Cov:   0.0784  -0.07826
              -0.07826  0.4225
    ===================================
    # This likelihood is from Cuceu, Farr, Lemos & Font-Ribera 2019
    # (1906.11628) and is based on data from Delubac et al. 2014 (1404.1801).
    # Written by James Farr (June 2019).
    
    BOSS_DR11_Lya_auto.data_directory      = data.path['data']
    BOSS_DR11_Lya_auto.file                = 'BOSS_DR11_scans/BOSS_DR11_scan_aux.txt'
    BOSS_DR11_Lya_auto.cf_scan             = 'BOSS_DR11_scans/BOSS_DR11_Lya_auto_scan.dat'
    BOSS_DR11_Lya_auto.rd_rescale          = 1.0
    
    BOSS_DR11_Lya_auto.transverse_fid      = 11.59
    BOSS_DR11_Lya_auto.parallel_fid        = 8.708
    '''

    def __init__(self):
        datafile = '/home/pub/ekli_work/cosmodata/BAO/BOSS_DR11_scans/BOSS_DR11_Lya_auto_scan.dat'
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
        interps = RectBivariateSpline(at,ap,grid,kx=1,ky=1)

        #Add the dictionary to the object.
        self.interp = interps
        
        return

    def loglikelihood(self, Hofz, DA, double rd=147.78):
        c = 2.99792458e5 #[km/s]
        zeff = 2.36
        transverse_fid = 11.59
        parallel_fid = 8.708

        #chisq = 0

        para = c/(Hofz(zeff)*rd) / parallel_fid
        if(np.isnan(para) or para<0.8 or para > 1.3):
            return -3e30

        trans = DA(zeff)/rd / transverse_fid
        if(np.isnan(trans) or trans <0.599999 or trans >1.5 ):
            return -3e30

        chisq = self.interp(trans, para)[0,0]
        
        return -0.5*chisq


