#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: test_test.py
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-05-02 19:38:34
#==================================

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

from gp_likelihood import GaussianProcesses as GP
from gp_likelihood import min_GaussianProcesses as mGP
from gp_likelihood import rec_GaussianProcesses as rGP
from gp_likelihood import rec_extern_GaussianProcesses as rexGP
from gapp import dgp

import time
from multiprocessing.dummy import Pool as ThreadPool

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#rec_dat = np.loadtxt('/home/ekli/myworks/cosmodata/Growth_tableII.txt', unpack=True)
rec_dat = np.loadtxt('/home/ekli/myworks/cosmodata/CC_31.txt', unpack=True)
X = rec_dat[0]
Y = rec_dat[1]
err = rec_dat[2]
cov = np.eye(len(err))*err**2

gpY = mGP(X,Y,cov)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
res = gpY.min_loglikelihood()
print(res)
print(res.x)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# reconst \int_0^x 1/mu(x) dx

sigf, lenf = res.x
gp = rexGP(X,Y,cov,sigf,lenf)

z = np.linspace(0,2,41)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
print('The start time is: %s'%time.asctime())
s_time = time.perf_counter()
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

Afunc = lambda x: 1
Bfunc = lambda x: 1
bet = 0.0
alp = -1.0

for zi in z:
    ff = gp.int_func(zi, Afunc, alp)
    #ff0 = gp.rec_mu_x(zi)**alp
    func = lambda x: gp.int_funcx(x, zi, Bfunc, alp)
    ff0, err = integrate.quad(func, 0.0, zi)
    #ff1 = gp.integrate_f_pq_x(Afunc, Bfunc, zi, 0.0, alp, bet)
    ff1, err = integrate.quad(gp.int_func, 0.0, zi, args=(Bfunc, alp))
    ff2 = gp.integrate_f_mu_x(zi,0.0, alp, bet)
    ff3 = gp.Cov_integrate_f_mu_x(zi, zi, 0.0, alp, bet)
    ff4 = gp.Cov_integrate_f_pq_x(Afunc,Bfunc,zi,zi,0.0,alp,bet)
    #ff4 = gp.int_funcx(zi, zi, Bfunc, alp)
    ff5 = gp.int_funcxy(zi, zi, Bfunc, alp)
    print("%10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e"
            %(ff, ff0, ff1, ff2, ff3, ff4, ff5))

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
e_time = time.perf_counter()
print('Using time %s s'%(e_time - s_time))
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
