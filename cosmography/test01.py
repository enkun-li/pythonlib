#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: test.py
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-04-18 00:08:40
#==================================

import numpy as np
import matplotlib.pyplot as plt
from cosmography import cosmolgoy_bk as cbk
from darksector import darksector as dks

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lcdm = cbk(H0=70,omm=0.3,omb=0.0462,MOD='lcdm', feedback=1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
edecdm = cbk(H0=67.27,omm=0.31,w0=-1.01,nps=np.array([0.01,0.01]),
        MOD='edecdm', feedback=1)

wx = lambda z: np.array([edecdm.EoS(zi) for zi in z])
omx = lambda z: np.array([edecdm.Omega_x(zi) for zi in z])

xx = np.linspace(-7.6,0,100)
zz = np.exp(-xx)-1

fig = plt.figure()
axs = fig.add_subplot(111)

axs.plot(1+zz, omx(zz), '--k')
axs.set_xscale('log')
axs.set_xlabel('1+z')
axs.set_ylabel('Omx(z)')

ax2 = axs.twinx()
ax2.plot(1+zz, wx(zz), '--r')
ax2.axhline(0.0, ls=':', color='k')
ax2.set_ylabel('w(z)')

#plt.show()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nps = np.array([-0.9,0.01])
IDE01 = cbk(H0=70,omm=0.3,omb=0.0462,nps=nps,MOD='IDE01', feedback=1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# nps 的几个参数分别对应 w_0, w_a, \xi_c, \xi_x
nps = np.array([-1.01,0.2, 0.0, 0.0])
IDE02 = cbk(H0=70,omm=0.3,omb=0.0462,nps=nps,MOD='IDE02', feedback=1)
#dk = dks(omegac=0.26,omegax=0.7,nps=nps,mod=5)
#dk.initial_rho_DM_DE()
#xx = dk.get_rho_DM_DE()
#xx = dk.w_x(0.3)

#print(xx)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
x = np.linspace(-7.6,0,50)
zz = np.exp(-x) -1

plt.figure()

hz = np.array([lcdm.Hofz(z) for z in zz])
plt.plot(zz+1, hz, label='LCDM')

hz = np.array([IDE01.Hofz(z) for z in zz])
plt.plot(zz+1, hz, label='IDE01')

hz = np.array([IDE02.Hofz(z) for z in zz])
plt.plot(zz+1, hz, label='IDE02')

plt.xscale('log')
plt.yscale('log')
plt.legend(loc='best')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plt.figure()

da = np.array([lcdm.angulardiameterdistance(z) for z in zz])
plt.plot(1+zz, da, label='LCDM')

da = np.array([IDE01.angulardiameterdistance(z) for z in zz])
plt.plot(1+zz, da, label='IDE01')

da = np.array([IDE02.angulardiameterdistance(z) for z in zz])
plt.plot(1+zz, da, label='IDE02')

plt.xscale('log')
plt.legend(loc='best')
plt.show()
