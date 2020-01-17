#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: likelihood.pyx
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-02-26 15:09:09
#==================================

cimport cython
import numpy as np
import array

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

#********************************
# loglikelihood for H0prior
#
cpdef double H0prior(Hth, Hobs, sigs):
    '''
    H0prior(Hth,Hobs,sigs)
    =========================
    param Hth: theortical H0
    param Hobs: observational data
    param sigs: observational sigma
    '''
    cdef double chisq
    chisq = -0.5*(Hth-Hobs)**2/sigs**2
    return chisq

#********************************
# loglikelihood for CC data
#
cdef class CC_loglikelihood:
    '''
    loglikelihood for CC data:
    \chi^2 = \sum_i (H_{theroy}(zi)-H_{obs}(zi))^2/sigma_i^2
    ============================
    Initialize likelihood:
    - Input: datafile
    -- datafile: the CC data filename
    '''
    cdef double[:] zcc, Hcc, ecc
    cdef int ncc
    
    def __cinit__(self, datafile):
        dats = np.loadtxt(datafile, unpack='True')
        self.ncc = np.shape(dats)[1]
        self.zcc = dats[0]
        self.Hcc = dats[1]
        self.ecc = dats[2]
        return

    cpdef double loglikelihood(self, Hofz):
        '''
        loglikelihood(Hofz)
        =======================
        Hofz: the hubble parameter function
        '''
        cdef double chi, chisq, Hth
        chisq = 0.0
        for i from 1 <= i < self.ncc:
            Hth = Hofz(self.zcc[i])
            if(np.isnan(Hth) or Hth <= 0):
                chisq = -1.0e30
                return chisq
            chi = (self.Hcc[i] - Hth)**2/self.ecc[i]**2
            chi += log(2*np.pi*self.ecc[i]**2)
            if(chi > 2.0e30):
                chisq = -1.0e30
                return chisq
            else:
                chisq += chi
        chisq = -chisq/2.0
        return chisq


#********************************
# loglikelihood for Pantheon data
#
cdef class SN_loglikelihood_dlum:
    '''
    loglikelihood for the Pantheon data:
    - chisq = (mu_o - mu_t)^T [C]^{-1} (mu_o -mu_t)
    - mu_t = 5 log10[ d_L ] + mu0
    - mu_o = m - M'
    - C_{ij} = \delta_{ij} sigma^2_{SN, i} + C^{sys}_{ij}
    -----------------------------------------------
    chisq = A - B^2/C + ln(C/2/pi)
    ===============================================
    Initialize likelihood:
    - Input: datafile, syscovfile
    -- datafile: the SN data filename [3 cloumns]
    -- syscovfile: the covariance matrix of systemtic [1 cloumns]
    '''
    cdef double[:] zsn, musn, esn
    cdef double[:,:] covsn, invcov
    cdef int nsn

    def __cinit__(self):
        cdef int i
        sndat = np.load('/home/ekli/myworks/cosmodata/lcparam_sn.npy')
        sncov = np.load('/home/ekli/myworks/cosmodata/syscov_sn.npy')
        self.nsn = np.shape(sndat)[1]
        self.zsn = sndat[0]
        self.musn = sndat[1]
        self.esn = sndat[2]
        self.covsn = sncov
        for i from 0<= i< self.nsn:
            self.covsn[i,i] += self.esn[i]**2
        self.invcov = np.linalg.inv(self.covsn)
        return

    cpdef double loglikelihood(self, dlum):
        '''
        loglikelihood(dlum)
        =========================
        dlum: the luminosity distance [d_L in [Mpc]]
        '''
        cdef int i
        cdef double zi, dl
        cdef double[:] difmu, Acov
        cdef double chisq, magA, magB, magC

        chisq = 0.0
        difmu = np.zeros(self.nsn)
        Acov = np.zeros(self.nsn)
        for i from 0<= i < self.nsn:
            zi = self.zsn[i]
            dl = dlum(zi)
            if(dl <= 0.0):
                chisq = -1.0e30
                return chisq
            else:
                difmu[i] = 5.0*log10(dl)-self.musn[i]
        Acov = np.matmul(difmu, self.invcov)
        magA = np.dot(Acov, difmu)
        magB = np.sum(Acov)
        magC = np.sum(self.invcov)
        chisq = magA - magB**2/magC + log(magC/np.pi)
        chisq = -chisq/2.0
        return chisq

    cpdef double argM(self, dlum):
        '''
        argM(dlum)
        =========================
        dlum: the luminosity distance [d_L in [Mpc]]
        '''
        cdef int i
        cdef double zi, dl
        cdef double[:] difmu, Acov
        cdef double chisq, magA, magB, magC

        chisq = 0.0
        difmu = np.zeros(self.nsn)
        Acov = np.zeros(self.nsn)
        for i from 0<= i < self.nsn:
            zi = self.zsn[i]
            dl = dlum(zi)
            if(dl <= 0.0):
                chisq = -1.0e30
                return chisq
            else:
                difmu[i] = 5.0*log10(dl)-self.musn[i]
        Acov = np.matmul(difmu, self.invcov)
        magA = np.dot(Acov, difmu)
        magB = np.sum(Acov)
        magC = np.sum(self.invcov)
        chisq = magB/magC
        return chisq

cdef class SN_loglikelihood:
    '''
    loglikelihood for the Pantheon data:
    - chisq = (mu_o - mu_t)^T [C]^{-1} (mu_o -mu_t)
    - mu_t = 5 log10[ d_L ] + mu0
    - mu_o = m - M'
    - C_{ij} = \delta_{ij} sigma^2_{SN, i} + C^{sys}_{ij}
    -----------------------------------------------
    chisq = A - B^2/C + ln(C/2/pi)
    ===============================================
    Initialize likelihood:
    - Input: datafile, syscovfile
    -- datafile: the SN data filename [3 cloumns]
    -- syscovfile: the covariance matrix of systemtic [1 cloumns]
    '''
    cdef double[:] zsn, musn, esn
    cdef double[:,:] covsn, invcov
    cdef int nsn

    def __cinit__(self):
        cdef int i
        sndat = np.load('/home/ekli/myworks/cosmodata/lcparam_sn.npy')
        sncov = np.load('/home/ekli/myworks/cosmodata/syscov_sn.npy')
        self.nsn = np.shape(sndat)[1]
        self.zsn = sndat[0]
        self.musn = sndat[1]
        self.esn = sndat[2]
        self.covsn = sncov
        for i from 0<= i< self.nsn:
            self.covsn[i,i] += self.esn[i]**2
        self.invcov = np.linalg.inv(self.covsn)
        return

    cpdef double loglikelihood(self, comdm):
        '''
        loglikelihood(comdm)
        =========================
        comdm: the comoving angular diameter distance [d_M in [Mpc]]
        d_L = (1+z) d_M
        '''
        cdef int i
        cdef double zi, dm
        cdef double[:] difmu, Acov
        cdef double chisq, magA, magB, magC

        chisq = 0.0
        difmu = np.zeros(self.nsn)
        Acov = np.zeros(self.nsn)
        for i from 0<= i < self.nsn:
            zi = self.zsn[i]
            dm = comdm(zi)
            if(dm <= 0.0):
                chisq = -1.0e30
                return chisq
            else:
                difmu[i] = 5.0*log10((1+zi)*dm)-self.musn[i]
        Acov = np.matmul(difmu, self.invcov)
        magA = np.dot(Acov, difmu)
        magB = np.sum(Acov)
        magC = np.sum(self.invcov)
        chisq = magA - magB**2/magC + log(magC/np.pi)
        chisq = -chisq/2.0
        return chisq
    
    cpdef double argM(self, comdm):
        '''
        argM(comdm)
        =========================
        comdm: the comoving angular diameter distance [d_M in [Mpc]]
        d_L = (1+z) d_M
        '''
        cdef int i
        cdef double zi, dm
        cdef double[:] difmu, Acov
        cdef double chisq, magA, magB, magC

        chisq = 0.0
        difmu = np.zeros(self.nsn)
        Acov = np.zeros(self.nsn)
        for i from 0<= i < self.nsn:
            zi = self.zsn[i]
            dm = comdm(zi)
            if(dm <= 0.0):
                chisq = -1.0e30
                return chisq
            else:
                difmu[i] = 5.0*log10((1+zi)*dm)-self.musn[i]
        Acov = np.matmul(difmu, self.invcov)
        magA = np.dot(Acov, difmu)
        magB = np.sum(Acov)
        magC = np.sum(self.invcov)        
        chisq = magB/magC
        return chisq

#********************************
# loglikelihood for DR12BAO
#
cdef class DR12BAO_loglikelihood:
    '''
    #BAO-only consensus results, Alam et al. 2016
    #https://arxiv.org/abs/1607.03155
    name = DR12BAO
    num_bao=6
    bao_measurements_file = sdss_DR12Consensus_bao.dat
    bao_measurements_file_has_error = F
    #stored values are actually scaled by r_fid, so include that here
    #rs_rescale = 1/147.78
    rs_rescale = .6766815537e-2
    bao_cov_file = BAO_consensus_covtot_dM_Hz.txt
    =====================================================
    * sdss_DR12Consensus_bao.dat:    
    -- 0.38 1512.39  DM_over_rs
    -- 0.38 81.2087  bao_Hz_rs
    -- 0.51 1975.22  DM_over_rs
    -- 0.51 90.9029  bao_Hz_rs
    -- 0.61 2306.68  DM_over_rs
    -- 0.61 98.9647  bao_Hz_rs
    =====================================================
    Input: path
    * path: the dir name where the cov file store
    '''
    cdef int nbao
    cdef double rs_rescale
    cdef double[:,:] baocov, invcov
    cdef int[:] type_bao
    cdef double[:] zbao, dbao

    def __cinit__(self, path):
        bao_cov_file = 'BAO_consensus_covtot_dM_Hz.txt'
        self.baocov = np.loadtxt(path+bao_cov_file, unpack=True)
        self.invcov = np.linalg.inv(self.baocov)
        self.nbao = 6
        self.rs_rescale = 0.6766815537e-2
        self.type_bao = array.array('i',[1,2,1,2,1,2])
        self.zbao = np.array([0.38, 0.38, 0.51, 0.51, 0.61, 0.61], \
                dtype=np.float64)
        self.dbao = np.array([1512.39, 81.2087, 1975.22, 90.9029, 2306.68, 98.9647], \
                dtype=np.float64)
        return

    cpdef double loglikelihood(self, dmz, Hofz, double rs=147.78):
        '''
        loglikelihood(dmz, Hofz, rs)
        ==============================
        The loglikelihood of DR12BAO
        - dmz:  the comoving angular diameter distance
        - Hofz: the Hubble parameter
        - rs:   the sound horizon at redshit z_drag [default: 147.78]
        '''
        cdef int i
        cdef double[:] difdbao
        cdef double chisq, dm, hz, rd

        rd = rs*self.rs_rescale

        difdbao = np.zeros(self.nbao)
        chisq = 0.0

        for i from 0<= i< self.nbao:
            dm = dmz(self.zbao[i])
            hz = Hofz(self.zbao[i])
            if(np.isnan(hz) or hz <=0.0  or dm <= 0):
                chisq = -1.0e30
                return chisq
            if(self.type_bao[i] == 1):
                difdbao[i] = dm/rd - self.dbao[i]
            elif(self.type_bao[i] == 2):
                difdbao[i] = hz*rd - self.dbao[i]
            else:
                print("No such type @DR12BAO")
                return chisq
        chisq = np.dot(np.matmul(difdbao, self.invcov), difdbao)
        chisq = -0.5*chisq
        return chisq

#******************************
# loglikelihood for DR14quasar
#
cdef class DR14quasar_loglikelihood:
    '''
    #DR14 quasar, https://arxiv.org/abs/1705.06373
    # Note alternative in http://arxiv.org/abs/1801.02656
    measurement_type = DV_over_rs
    zeff = 1.52
    bao_measurement = 26.086 1.150
    ==============================
    no input
    '''
    cdef double zeff, dvrs, err
    
    def __cinit__(self):
        self.zeff = 1.52
        self.dvrs = 26.086
        self.err = 1.150
        return

    cpdef double loglikelihood(self, bao_dv, double rs=147.78):
        '''
        loglikelihood(bao_dv, rs) = DV_over_rs
        ======================================
        The loglikelihood of sdss DR14 qiasar bao:
        - bao_dv: the effective distance of bao volume
        - rs: the sound speed at redshit z_drag
        '''
        cdef double chisq, dv
        dv = bao_dv(self.zeff)
        if(np.isnan(dv)):
            chisq = -1.0e30
            return chisq
        chisq = (dv/rs - self.dvrs)**2/self.err**2
        chisq = -0.5*chisq        
        return chisq

#******************************
# loglikelihood for 6DF
#
cdef class DF6_loglikelihood:
    '''
    # http://arxiv.org/abs/1106.3366
    name = 6DF
    measurement_type = rs_over_DV
    #rs_rescale = 153.9/149.8 
    rs_rescale= 1.027369826
    zeff= 0.106
    bao_measurement = 0.336 0.015
    ==============================
    no input
    '''
    cdef double zeff, dvrs, err, rs_rescale
    
    def __cinit__(self):
        self.zeff = 0.106
        self.dvrs = 0.336
        self.err = 0.015
        self.rs_rescale = 1.027369826
        return

    cpdef double loglikelihood(self, bao_dv, double rs=147.78):
        '''
        loglikelihood(bao_dv, rs) = rs_over_DV
        ======================================
        The loglikelihood of sdss 6DF data:
        - bao_dv: the effective distance of bao volume
        - rs: the sound speed at redshit z_drag
        '''
        cdef double chisq, dv, rd
        rd = rs*self.rs_rescale
        dv = bao_dv(self.zeff)
        if(np.isnan(dv)):
            chisq = -1.0e30
            return chisq
        chisq = (rd/dv - self.dvrs)**2/self.err**2
        chisq = -0.5*chisq        
        return chisq

#**********************************************
# distance priors from Planck Final Release
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cdef class distance_loglikelihood:
    '''
    The distance priors form the finally released Planck TT,TE,EE+lowE data in
    2018.
    =============================
    r = 1.750235
    la = 301.4707
    omegabh2 = 0.02235976
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    distance invcov:
    # PLA 2018 base_plikHM_TTTEEE_lowE
    # R l_a a% omegabh2
      94392.3971 -1360.4913   1664517.2916
      -1360.4913   161.4349      3671.6180
    1664517.2916  3671.6180  79719182.5162
    '''
    cdef double R, la, ombh2
    cdef double[:,:] invC
    cdef double c

    def __cinit__(self):
        self.R = 1.750235
        self.la = 301.4707
        self.ombh2 = 0.02235976
        invC = np.array([
            [  94392.3971, -1360.4913,  1664517.2916],
            [  -1360.4913,   161.4349,     3671.6180],
            [1664517.2916,  3671.6180, 79719182.5162]])
        self.c = 2.99792458e5 #[km/s]
        return
    
    cpdef double loglikelihood(self, DA, double theta, double omb, double omm, double h0):
        '''
        loglikelihood(DA, theta, omb, omm, h0)
        ==========================================
        DA(z): angular diameter distance function
        theta: rs/DA at zstar
        omb: Omega_b - energy density of baryons
        omm: Omega_m - energy density of matter
        h0: Hubble constant
        '''
        cdef double g1, g2, zstar, R, la
        cdef double ombh2, ommh2
        cdef double[:] dx
        cdef double chisq

        ombh2 = omb*(h0/100.0)**2
        ommh2 = omm*(h0/100.0)**2
        g1 = 0.0783 *ombh2**(-0.238)/(1.0 +39.5 *ombh2**0.763)
        g2 = 0.560/(1.0 +21.1*ombh2**1.81)
        zstar = 1048.0*(1.0 +0.00124 *ombh2**(-0.738))*(1.0+g1*ommh2**g2)
        
        R = h0*omm**0.5*DA(zstar)*(1.0+ zstar)/self.c
        la = np.pi/theta

        dx = np.array([self.R-R, self.la-la, self.ombh2-ombh2])
        chisq = -0.5*np.dot(dx, np.matmul(self.invC, dx))
        return chisq

#**********************************************
# distance priors from PLA release
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# la Ra and z* all could be derive from 
# zstar, rstar, 100thetastar, DAstar
#
#cdef cl
