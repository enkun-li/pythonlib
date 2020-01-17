#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: test_ext.py
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-05-01 16:44:08
#==================================

import numpy as np
import matplotlib.pyplot as plt
from gp_likelihood import GaussianProcesses as GP
from gp_likelihood import min_GaussianProcesses as mGP
from gp_likelihood import rec_GaussianProcesses as rGP
from gp_likelihood import rec_extern_GaussianProcesses as rexGP
from gapp import dgp

import time
from multiprocessing import Pool, cpu_count

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
def get_rec_int_dzdmu(z,X,Y,cov,sigf,lenf):
    print('The start time is: %s'%time.asctime())
    s_time = time.perf_counter()

    gp = rexGP(X,Y,cov,sigf,lenf)
    
    rec = np.array([gp.integrate_one_over_mu(zi,0.0) for zi in z])
    sigsq = np.array([gp.Cov_integrate_one_over_mu(zi,zi,0) for zi in z])
    
    e_time = time.perf_counter()
    print('Using time %s s'%(e_time - s_time))
    return rec, sigsq

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
def get_integrate_f_mu_x(z,X,Y,cov,sigf,lenf):
    print('The start time is: %s'%time.asctime())
    s_time = time.perf_counter()
    
    gp = rexGP(X,Y,cov,sigf,lenf)
    
    rec = np.array([gp.integrate_f_mu_x(zi,0.0,-1,0) for zi in z])
    sigsq = np.array([gp.Cov_integrate_f_mu_x(zi,zi,0,-1,0) for zi in z])
    
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    e_time = time.perf_counter()
    print('Using time %s s'%(e_time - s_time))

    return rec, sigsq

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
def get_integrate_f_pq_x(z,X,Y,cov,sigf,lenf):
    print('The start time is: %s'%time.asctime())
    s_time = time.perf_counter()

    gp = rexGP(X,Y,cov,sigf,lenf)

    Afunc = lambda x: 1.0
    Bfunc = lambda x: 1.0
    bet = 0.0
    alp = -1.0
    
    rec = np.array([gp.integrate_f_pq_x(Afunc,Bfunc,zi,0.0,alp,bet) for zi in z])
    sigsq = np.array([gp.Cov_integrate_f_pq_x(Afunc,Bfunc,zi,zi,0.0,alp,bet) for zi in z])
    
    e_time = time.perf_counter()
    print('Using time %s s'%(e_time - s_time))
    return rec, sigsq

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
def output_fig(rec,sigsq,fig=1, sty='--', color='red'):
    ff = rec*300
    sig = np.sqrt(sigsq)*300
    
    plt.figure(fig)
    
    plt.plot(z, ff, ls='--', color=color)
    plt.fill_between(z, ff+sig, ff-sig, color=color, alpha=0.5)
    return

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
def multi_get_integrate_f_pq_x(z,X,Y,cov,sigf,lenf):
    print('The start time is: %s'%time.asctime())
    s_time = time.perf_counter()

    gp = rexGP(X,Y,cov,sigf,lenf)

    Afunc = lambda x: 1.0
    Bfunc = lambda x: 1.0
    bet = 0.0
    alp = -1.0

    n = z.shape[0]

    #~~~~~~~~~~~~~~~~
    cores = cpu_count()
    pool1 = Pool(cores)
    pool2 = Pool(cores)
    p1_list = []
    p2_list = []

    for i in range(n):
        pp = pool1.apply_async(gp.integrate_f_pq_x,
                args=(Afunc,Bfunc,z[i],0.0,alp,bet) )
        ps = pool2.apply_async(gp.Cov_integrate_f_pq_x,
                args=(Afunc,Bfunc,z[i],z[i],0.0,alp,bet) )
        p1_list.append(pp)
        p2_list.append(ps)
    #result1 = np.array([pp.get() for pp in p1_list])
    #result2 = np.array([ps.get() for ps in p2_list])

    rec = [pp.get() for pp in p1_list]
    sigsq = p2_list

    pool1.close()
    pool2.close()
    pool1.join()
    pool2.join()
    
    e_time = time.perf_counter()
    print('Using time %s s'%(e_time - s_time))
    return rec, sigsq


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__=='__main__':
    # reconst \int_0^x 1/mu(x) dx
    datfile = '/home/ekli/myworks/cosmodata/CC_31.txt'
    X,Y,cov = readin_dat(datfile)
    sigf, lenf = get_hyper_pars(X,Y,cov)
    
    z = np.linspace(0,2,41)
    rec0, sigsq0 = get_rec_int_dzdmu(z,X,Y,cov,sigf,lenf)
    rec1, sigsq1 = get_integrate_f_mu_x(z,X,Y,cov,sigf,lenf)
    rec2, sigsq2 = get_integrate_f_pq_x(z,X,Y,cov,sigf,lenf)
    #rec3, sigsq3 = multi_get_integrate_f_pq_x(z,X,Y,cov,sigf,lenf)

    output_fig(rec0,sigsq0,fig=1,sty='--',color='red')
    output_fig(rec1,sigsq1,fig=2,sty='--',color='blue')
    output_fig(rec2,sigsq2,fig=3,sty='--',color='green')
    #output_fig(rec3,sigsq3,fig=4,sty='--',color='grey')

    plt.show()
