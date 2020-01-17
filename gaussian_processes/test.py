#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: test.py
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-04-29 10:33:21
#==================================

import numpy as np
import matplotlib.pyplot as plt
from gp_likelihood import GaussianProcesses as GP
from gp_likelihood import min_GaussianProcesses as mGP
from gp_likelihood import rec_GaussianProcesses as rGP
from gp_likelihood import rec_extern_GaussianProcesses as rexGP
from gapp import dgp

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def gapp_Y(X, Y, err):
    # Gaussian process
    zmin = np.min(X)
    zmax = np.max(X)
    g = dgp.DGaussianProcess(X,Y,err, cXstar=(zmin,zmax,40),theta=[0.5,0.5])
    (rec,theta) = g.gp()
    (drec,theta) = g.dgp(thetatrain='False')
    (d2rec,theta) = g.d2gp()    
    
    #plt.show()
    return theta, rec, drec, d2rec
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#rec_dat = np.loadtxt('/home/ekli/myworks/cosmodata/Growth_tableII.txt', unpack=True)
data_path = '/home/ekli/myworks/cosmodata/'
rec_dat = np.loadtxt(data_path+'OHD_51.txt', unpack=True)
X = rec_dat[0]
Y = rec_dat[1]
err = rec_dat[2]
cov = np.eye(len(err))*err**2

gpY = mGP(X,Y,cov)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# using different bonds
#
#Ymax = np.max(Y)
#zmax = np.max(X)
#
##res = gpY.min_loglikelihood(2*Ymax, 2*zmax)
##print(res.x)
#
#for i in range(10):
#    res = gpY.min_loglikelihood(2*(i+1)*Ymax,2*(i+1)*zmax)
#    print('%3s-th is: '%(i+1), res.x)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
res = gpY.min_loglikelihood()
print(res)
print(res.x)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sigf, lenf = res.x
gp = rGP(X,Y,cov,sigf,lenf)

z = np.linspace(0,2,41)

recfx = gp.rec_mu(z)
covfx = gp.rec_covarianve(z)
sig = np.array([np.sqrt(covfx[i,i]) for i in range(len(z))])

recdfx = gp.rec_dmu(z)
covdfx = gp.rec_dfdfcovarianve(z)
sigdf = np.array([np.sqrt(covdfx[i,i]) for i in range(len(z))])

#recddfx = gp.rec_ddmu(z)
#covddfx = gp.rec_ddfddfcovarianve(z)
#sigddf = np.array([np.sqrt(covddfx[i,i]) for i in range(len(z))])

recddfx = gp.rec_ddmu(z)
covddfx = gp.rec_diff_covarianve(z,dx=2,dy=2)
sigddf = np.array([np.sqrt(covddfx[i,i]) for i in range(len(z))])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#sigf, lenf = res.x
#xx = -np.log(1+X)
#gp = rGP(xx,Y,cov,sigf,lenf)
#
#z = np.linspace(0,2,41)
#x = -np.log(1+z)
#
#recfx = gp.rec_mu(x)
#covfx = gp.rec_covarianve(x)
#
#sig = np.array([np.sqrt(covfx[i,i]) for i in range(len(z))])
#
#plt.plot(z, recfx, '--r')
#plt.errorbar(z, recfx, sig, color='red')
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
theta, rec, drec, d2rec = gapp_Y(X,Y,err)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plt.figure()

#-----------------------------------------
plt.subplot(3,1,1)

plt.plot(z, recfx, '--r')
#plt.errorbar(z, recfx, sig, color='red')
plt.fill_between(z, recfx+sig, recfx-sig, color='red', alpha=0.5)

plt.errorbar(X, Y, yerr=err, 
        color='grey', fmt='o', markersize=8, lw=2, capthick=1,
        capsize=1+2*2,elinewidth=2,label='RSD Data Points',
        markerfacecolor='white')

plt.plot(rec[:,0], rec[:,1], ':b')
plt.errorbar(rec[:,0], rec[:,1], rec[:,2], color='blue')

plt.xlabel('z')
plt.ylabel(r'$f\sigma_8(z)$')

#-----------------------------------------
plt.subplot(3,1,2)

plt.plot(z, recdfx, '--r')
plt.fill_between(z, recdfx+sigdf, recdfx-sigdf, color='red', alpha=0.5)

plt.plot(drec[:,0], drec[:,1], ':b')
plt.errorbar(drec[:,0], drec[:,1], drec[:,2], color='blue')

plt.xlabel('z')
plt.ylabel(r'$df\sigma_8(z)/dz$')

#------------------------------------------
plt.subplot(3,1,3)

plt.plot(z, recddfx, '--r')
plt.fill_between(z, recddfx+sigddf, recddfx-sigddf, color='red', alpha=0.5)

plt.plot(d2rec[:,0], d2rec[:,1], ':b')
plt.errorbar(d2rec[:,0], d2rec[:,1], d2rec[:,2], color='blue')

plt.xlabel('z')
plt.ylabel(r'$d2f\sigma_8(z)dz2$')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# reconst \int_0^x 1/mu(x) dx

sigf, lenf = res.x
gp = rexGP(X,Y,cov,sigf,lenf)

z = np.linspace(0,2,41)

rec_int_dzdmu = np.array([gp.integrate_one_over_mu(zi,0.0) for zi in z])
sig_int_dzdmu = np.array([gp.Cov_integrate_one_over_mu(zi,zi,0) for zi in z])

ff = rec_int_dzdmu
sig = np.sqrt(sig_int_dzdmu)

plt.figure()

plt.plot(z, ff, '--r')
plt.fill_between(z, ff+sig, ff-sig, color='green', alpha=0.5)


plt.show()
