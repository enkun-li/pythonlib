#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: test.py
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-01-13 22:49:17
#==================================

import numpy as np
import mcmc as mc
import matplotlib.pyplot as plt
import scipy.optimize as op

from likelihood import sn_likelihood as snlike

path = '/home/ekli/myworks/Model_indepent/dat2fig/'

snl = snlike(path+'outputdat/lcparam_sn.npy',\
        path+'outputdat/syscov_sn.npy',\
        path+'outputdat/rec_dc_sn.npy',\
        path+'outputdat/rec_covdcdc_sn.npy')

#snl.setup_pars(70.0, 0.0)

amag = snl.amag(70,0.0)

print(amag[0], amag[1], amag[2])
print('M is : %s'%(amag[1]/amag[2]))

for i in range(1):
    H0 = np.random.uniform(60.0,80.0)
    omk = np.random.uniform(-1,1)
    print('='*60)
    print('H0 is: %10.3f'%H0)
    print('omk is: %10.3f'%omk)
    print('lnlike is: %10.3f'%snl.lnlike_sn(H0, omk))

print('Use the likelihood to find the best fit parameters')
def lnprob(theta):
    H0, omk = theta
    return snl.lnlike_sn(H0,omk)

nll = lambda *args: -lnprob(*args)

iniH0 = 70.0
iniomk = 0.0
result = op.minimize(nll, [iniH0, iniomk])

print(result["x"])
