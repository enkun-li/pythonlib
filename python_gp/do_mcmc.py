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
import emcee
import time

from likelihood import sn_likelihood as snlike

path = '/home/ekli/myworks/Model_indepent/dat2fig/'

snl = snlike(path+'outputdat/lcparam_sn.npy',\
        path+'outputdat/syscov_sn.npy',\
        path+'outputdat/rec_dc_sn.npy',\
        path+'outputdat/rec_covdcdc_sn.npy')

print('Use the likelihood to find the best fit parameters')
#def lnprior(theta):
#    H0, omk = theta
#    if(60<H0<80 and -1<omk<1):
#        return 0.0
#    return -np.inf

#def lnprob(theta):
#    H0, omk = theta
#    lp = lnprior(theta)
#    if(not np.isfinite(lp)):
#        return -np.inf
#    return lp+snl.lnlike_sn(H0,omk)

initheta = [70.0, 0.0]
theta = np.zeros(2)

#for i in range(10):
#    theta[0] = 65 + (i-1)/9.0*(80-65.0)
#    for j in range(10):
#        theta[1] = -1+(i-1)/9.0*2.0
#        print('H0 = %f5.2, omk = %5.2f, like = %20.8f'\
#                %(theta[0], theta[1], snl.lnprob(theta)))

pos = initheta + 1e-4*np.random.randn(100, 2)
nwalkers, ndim = pos.shape
nsteps = 3000

sampler = emcee.EnsembleSampler(nwalkers, ndim, snl.lnprob)

print('do the mcmc')
start = time.time()
sampler.run_mcmc(pos, nsteps, progress=True)
end = time.time()
print('Serial took {0:.1f} seconds'.format(end -start))

print('Save the sample')
#np.save('sampler', sampler)
samples = sampler.chain[:,50:,:].reshape((-1,ndim))
np.save('samples', samples)
np.save('samples_%s'%nsteps, samples)
