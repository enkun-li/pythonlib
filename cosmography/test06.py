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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# nps 的几个参数分别对应 w_0, w_a, \xi
plt.figure(1)
plt.subplot(2,1,1)
x = np.linspace(-7.6,0,50)
zz = np.exp(-x) -1
hz = np.array([lcdm.Hofz(z) for z in zz])
plt.plot(zz+1, hz, label='LCDM')

plt.subplot(2,1,2)
zs = np.linspace(0,2,50)
hz = np.array([lcdm.Hofz(z) for z in zs])
plt.plot(zs+1, hz, label='LCDM')

for i in range(10):
    xi = np.random.uniform(-0.1,0.1)
    nps3 = np.array([-1.1,0.0, xi])
    IDE03 = cbk(H0=70,omm=0.3,omb=0.0462,nps=nps3,MOD='IDE03', feedback=1)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
    plt.figure(1)
    plt.subplot(2,1,1)
    hz = np.array([IDE03.Hofz(z) for z in zz])
    plt.plot(zz+1, hz, label='IDE03- %5.2f'%xi)
    
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc='best')
    
    #~~~~~~~~~~~~~~~~
    plt.subplot(2,1,2)
    
    hz = np.array([IDE03.Hofz(z) for z in zz])
    plt.plot(zz+1, hz, label='IDE03- %5.2f'%xi)
    
    plt.legend(loc='best')

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#plt.figure()
#plt.subplot(2,1,1)
#
#x = np.linspace(-7.6,0,50)
#zz = np.exp(-x) -1
#
#da = np.array([lcdm.angulardiameterdistance(z) for z in zz])
#plt.plot(1+zz, da, label='LCDM')
#
#da = np.array([IDE03.angulardiameterdistance(z) for z in zz])
#plt.plot(1+zz, da, label='IDE03')
#
#plt.xscale('log')
#plt.legend(loc='best')
#
##~~~~~~~~~~~~~~~~~
#plt.subplot(2,1,2)
#zz = np.linspace(0,2,50)
#
#da = np.array([lcdm.angulardiameterdistance(z) for z in zz])
#plt.plot(1+zz, da, label='LCDM')
#
#da = np.array([IDE03.angulardiameterdistance(z) for z in zz])
#plt.plot(1+zz, da, label='IDE03')
#
#plt.xscale('log')
#plt.legend(loc='best')
plt.show()
