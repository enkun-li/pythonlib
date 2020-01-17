#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: test.py
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-06-21 19:45:04
#==================================

from my_gapp import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

#==========================================
def lcdm_hz(z,omm=0.3,h0=70):
    return h0*np.sqrt(omm*(1+z)**3+1-omm)

def lcdm_dz(z,omm=0.3,h0=70):
    func = lambda x: 1/lcdm_hz(x,omm,h0)
    dz, err = quad(func, 0.0, z)
    return h0*dz

def lcdm_mu_p(z,omm=0.3,h0=70):
    try:
        n = len(z)
    except:
        n = 1
    if (n==1):
        return 5*np.log10((1+z)*lcdm_dz(z,omm,h0))
    dz = np.array([lcdm_dz(zi,omm,h0) for zi in z])
    return 5*np.log10((1+z)*dz)

def M_prime(mudl, mobs, sig):
    B = np.sum((mudl-mobs)/sig**2)
    C = np.sum(1/sig**2)
    return -B/C 

#==========================================

sn_file = bytes("/home/ekli/myworks/cosmodata/sn_full_dat.txt", "utf8")
sn_cov = bytes("/home/ekli/myworks/cosmodata/sn_full_cov.txt", "utf8")

sndat = np.loadtxt("/home/ekli/myworks/cosmodata/sn_full_dat.txt", unpack=True)

#Mb = sndat[1] - lcdm_mu_p(sndat[0])
#mMb = np.mean(Mb)
mudl = lcdm_mu_p(sndat[0])
mMb = M_prime(mudl, sndat[1], sndat[2])
print(mMb)

#===========================================

N = 1048

initial_gapp(N, sn_file)
initial_gapp_cov(N, sn_cov)
#setup_gapp(39.21533362, 12.97267722) # for lz
setup_gapp(24.49128409, 0.14743146) # for z

lz = np.linspace(-5.,1.0, 51)
zz = np.exp(lz)

mu = np.array([rec_mu(z) for z in zz])
sig = np.array([np.sqrt(rec_covariance(z,z) ) for z in zz])

dz = 1/(1+zz)*np.exp(np.log(10)/5*(mu))
lc_dz = np.array([lcdm_dz(z) for z in zz])

#============================
plt.figure()

plt.errorbar(sndat[0], sndat[1], sndat[2], 
        fmt='.', markerfacecolor='white',
        color='red', alpha=0.1)

plt.plot(zz, mu, '--k')
plt.fill_between(zz, mu+sig, mu-sig, color='blue', alpha=0.5)

plt.plot(zz, lcdm_mu_p(zz)+mMb, '-g')

plt.xscale('log')

#=============================
plt.figure()
plt.plot(zz, dz/10**(0.2*mMb), '--')
plt.plot(zz, lc_dz, '-k')

#=============================
plt.figure()
plt.plot(zz, 5*np.log10(dz/lc_dz), '.')

plt.show()


