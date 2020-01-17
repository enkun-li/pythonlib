#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: test_sn.py
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-04-30 09:02:16
#==================================

import numpy as np
import matplotlib.pyplot as plt
from gp_likelihood import GaussianProcesses as GP
from gp_likelihood import min_GaussianProcesses as mGP
from gp_likelihood import rec_GaussianProcesses as rGP
from gapp import dgp, gp

#sndat = np.load('/home/ekli/myworks/cosmodata/lcparam_sn.npy')
sndat = np.loadtxt('/home/ekli/myworks/cosmodata/sn_full_lnz_dat.txt', unpack=True)
sncov = np.load('/home/ekli/myworks/cosmodata/syscov_sn.npy')

X = sndat[0]
Y = sndat[1]
sig = sndat[2]
Cov = sncov + np.eye(len(X))*sig**2
#Cov = np.eye(len(X))*sig**2

#gpsn = mGP(X,Y,Cov)

#res = gpsn.min_loglikelihood()
#print(res)

xmin = np.min(X)
xmax = np.max(X)
nstar = 51

initheta = [np.max(Y), 5]
g = gp.GaussianProcess(X,Y, Cov, cXstar=(xmin,xmax,nstar))

rec, theta = g.gp(theta=initheta)

print(theta)
np.save('sn_rec', rec)


