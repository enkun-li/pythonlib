#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: test_fs8.py
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-05-18 22:17:53
#==================================

import matplotlib.pyplot as plt
import numpy as np
import time
from gp_likelihood import min_GaussianProcesses as mGP
from general_fsigma import GaussianProcesses as fs8GP

def readin_dat(datfile):
    #rec_dat = np.loadtxt('/home/ekli/myworks/cosmodata/Growth_tableII.txt', unpack=True)
    rec_dat = np.loadtxt(datfile, unpack=True)
    X = rec_dat[0]
    Y = rec_dat[1]
    err = rec_dat[2]
    cov = np.eye(len(err))*err**2
    return X, Y, cov

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def get_hyper_pars(X,Y,cov):
        
    gpY = mGP(X,Y,cov)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    res = gpY.min_loglikelihood()
    print(res)
    print(res.x)
    return res.x

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
def get_rec_fsigma(zz,X,Y,cov,sigf,lenf):
    print('The start time is: %s'%time.asctime())
    s_time = time.perf_counter()

    gp = fs8GP(X,Y,cov,zz,sigf,lenf)

    gam = 6.0/11.0

    recint1 = gp.integrate_fx_over_mu(gam)
    recint2 = gp.integrate_fx_over_mu_covx(gam)
    recint3 = gp.integrate_dbfx_over_mu_covxy(gam)

    rec_cov = gp.cov_over_mumu_mat()
    rech = gp.rec_mu_arr()

    n = rech.shape[0]
    cov = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            cov[i,j] = rec_cov[i,j] - recint2[i,j] -recint2[j,i] + recint3[i,j]
    
    e_time = time.perf_counter()
    print('Using time %s s'%(e_time - s_time))
    return recint1, recint2, recint3, rech, rec_cov, cov

def get_fsigma8_theory(sig8,omm,h0,gam,zz,X,Y,cov,sigf,lenf):
    print('The start time is: %s'%time.asctime())
    s_time = time.perf_counter()

    gp = fs8GP(X,Y,cov,zz,sigf,lenf)

    gp.initial_theory_fs8(sig8,omm,h0,gam)

    fs8 = gp.return_fs8_theory()
    cov = gp.return_fs8cov_theory()

    e_time = time.perf_counter()
    print('Using time %s s'%(e_time - s_time))
    return fs8, cov

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def test_get_rec_fsigma(z,X,Y,cov,sigf,lenf):
    r1,r2,r3,rh, rc, cov = get_rec_fsigma(z,X,Y,cov,sigf,lenf)
    
    plt.figure()
    plt.plot(z, r1, '--')

    plt.figure()
    im = plt.imshow(r2)
    plt.colorbar(im)

    plt.figure()
    im = plt.imshow(r3)
    plt.colorbar(im)

    plt.figure()
    plt.plot(z,rh, '-r')

    plt.figure()
    im = plt.imshow(rc)
    plt.colorbar(im)
    
    plt.figure()
    im = plt.imshow(cov)
    plt.colorbar(im)

    plt.show()
    return

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def test_get_fsigma8_theory(sig8,omm,h0,gam,z,X,Y,cov,sigf,lenf):
    fs8, cov = get_fsigma8_theory(sig8,omm,h0,gam,z,X,Y,cov,sigf,lenf)

    plt.figure()
    plt.plot(z, fs8)

    plt.figure()
    im = plt.imshow(cov)
    plt.colorbar(im)

    plt.show()
    return

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def test_1px_over_mu(zz,X,Y,cov,sigf,lenf):
    print('The start time is: %s'%time.asctime())
    s_time = time.perf_counter()
    
    gp = fs8GP(X,Y,cov,zz,sigf,lenf)

    hz = gp.rec_mu_arr()
    
    e_time = time.perf_counter()
    print('Using time %s s'%(e_time - s_time))
    
    plt.figure()
    plt.plot(zz, hz)

    plt.figure()
    plt.plot(zz, 70**2*(1+zz)**3/hz/hz)
    plt.show()
    return

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def test_int_1px_over_mu(zz,X,Y,cov,sigf,lenf):
    print('The start time is: %s'%time.asctime())
    s_time = time.perf_counter()
    
    gp = fs8GP(X,Y,cov,zz,sigf,lenf)

    gam = 6.0/11
    
    hz = gp.rec_mu_arr()
    inthz = gp.integrate_fx_over_mu(gam)
    
    e_time = time.perf_counter()
    print('Using time %s s'%(e_time - s_time))
    
    plt.figure()
    plt.plot(zz, 1/np.exp(inthz) )
    
    plt.figure()
    plt.plot(zz, (1+z)**(3*gam)/hz/hz/np.exp(inthz) )

    plt.show()

    return

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__=='__main__':
    # reconst \int_0^x 1/mu(x) dx
    datfile = '/home/ekli/myworks/cosmodata/OHD_51.txt'
    X,Y,cov = readin_dat(datfile)
    sigf, lenf = get_hyper_pars(X,Y,cov)

    z = np.linspace(0,2,41)

    sig8 = 0.8
    omm = 0.3
    h0 = 70
    gam = 6.0/11

    test_get_fsigma8_theory(sig8,omm,h0,gam,z,X,Y,cov,sigf,lenf)
    #test_1px_over_mu(z,X,Y,cov,sigf,lenf)
    #test_int_1px_over_mu(z,X,Y,cov,sigf,lenf)
    
    


