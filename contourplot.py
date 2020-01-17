#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: contourplot.py
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-04-03 21:13:41
#==================================
'''
This file is used to generate CMB shift datas, corresponding matrix covariance
matrix from the already existed MCMC chains.
You can change it to include more parameters.
'''

import numpy as np
import matplotlib.pyplot as plt
from getdist import plots as gplt
from getdist import MCSamples, loadMCSamples
import sys

c_light = 2.99792458e8 #[m/s]
info = {'cmbruler': ['zstar', 'rstar', 'thetastar', 'DAstar', 'lA', 'R', 'omegabh2'],
        'derive': ['lA', 'R']
        }
#===================================
def InvMat(cov):
    num = cov.shape[0]
    L = np.linalg.cholesky(cov)
    B = np.linalg.solve(L, np.eye(num))
    invM = np.linalg.solve(np.transpose(L), B)
    return invM

#===================================
def triangleplot(samples, pinfo='cmbruler'):
    g = gplt.getSubplotPlotter(width_inch=9.0)
    g.settings.axes_fontsize = 12.0
    g.settings.figure_legend_frame = False
    g.settings.num_plot_contours = 2
    pars = info[pinfo]
    g.triangle_plot(samples, pars)
    gplt.plt.show()
    return

def outputdats(samples,pinfo='cmbruler'):
    roots = samples.getName()
    ff = open(roots+'_'+pinfo+'.txt', 'w')
    print(roots)
    ff.writelines(['#', roots, '\n\n'])
    print(info[pinfo])
    ff.write('num_dat = %s\n\n'%len(info[pinfo]))
    means = samples.getMeans()
    covs = samples.cov()
    corrs = samples.corr()
    pname = samples.getParamNames()

    conts = [pname.numberOfName(ss) for ss in info[pinfo]]

    cov = np.array([[covs[i,j] for i in conts] for j in conts])
    corr = np.array([[corrs[i,j] for i in conts] for j in conts])
    
    ff.writelines(['#Parameter'.ljust(12), 'best-fit'.center(18),
        'mean'.center(18), '68% limits'.center(18), '\n'])
    for ss, cont in zip(info[pinfo], conts):
        pp = pname.parWithName(ss)
        bf = pp.bestfit_sample
        mu = means[cont]
        sig = np.sqrt(covs[cont, cont])
        lab = pp.latexLabel()
        print(ss, bf, mu, sig, lab)
        ff.writelines([ss.ljust(12), '%18.8f'%bf, '%18.8f'%mu, '%18.8f'%sig,
            '    %s'%lab, '\n'])
    print(conts)
    
    print(cov)
    num = len(info[pinfo])
    ff.writelines(['\n', '#Covariance Matrix', '\n'])
    for i in range(num):
        ff.writelines(['%18.8e'%cov[i,j] for j in range(num)])
        ff.write('\n')
    
    print(corr)
    ff.writelines(['\n', '#Correlation Matrix', '\n'])
    for i in range(num):
        ff.writelines(['%18.8f'%corr[i,j] for j in range(num)])
        ff.write('\n')
        
    invC = InvMat(cov)
    ff.writelines(['\n', '#Inverse Covariance Matrix', '\n'])
    for i in range(num):
        ff.writelines(['%18.8e'%invC[i,j] for j in range(num)])
        ff.write('\n')

    ff.close()
    return


def outputcorrelation(roots,samples,pinfo='cmbruler'):
    #covmat = samples.getCovMat()
    #covmat.plot()
    #gplt.plt.show()
    cormat = samples.getCorrelationMatrix()
    #cm = plt.imshow(cormat)
    #plt.colorbar(cm)

    pname = samples.getParamNames()
    #means = samples.getMeans()

    labs = info[pinfo]
    num = len(labs)
    conts = []
    cont = 0
    ff = open(roots+'_'+pinfo+'.txt', 'w')
    ff.write('# '+roots+'\n\n')
    ff.write('%s %s %s\n'%(r'# Parameter'.ljust(12), r'best fit'.center(18), '68% limits'.center(18)))

    l1 = np.zeros((num, 2))
    i = 0

    for aa in pname.list():
        for lb in labs:
            if aa == lb:
                print('%2s %3s %s'%(cont, ':', aa))
                conts.append(cont)
                pp = pname.parWithName(aa)
                bf = pp.bestfit_sample
                lim1 = samples.getInlineLatex(param=aa, limit=1)
                l1s = lim1.split(sep='=')[-1]
                l1ss = l1s.split(sep='\\pm')
                l1[i] = [float(l1ss[0]), float(l1ss[1])]
                #print(l1s, l1[i])
                i += 1
                ff.write('%s = %16.8f %s\n'%(aa.ljust(10), bf, l1s))
        cont +=1
    print(conts)

    corrs = np.zeros((num,num))
    cov = np.zeros((num,num))
    cij = lambda rij,si,sj: rij*si*sj

    ff.write('\n# correlation matrix \n')
    for i in range(num):
        for j in range(num):
            corrs[i,j] = cormat[conts[i], conts[j]]
            cov[i,j] = cij(corrs[i,j],l1[i,1], l1[j,1])
            ff.write('%16.8f'%corrs[i,j])
        ff.write('\n')

    ff.write('\n# covariance matrix \n')
    for i in range(num):
        for j in range(num):
            ff.write('%16.8f'%cov[i,j])
        ff.write('\n')
    
    #invC = np.linalg.inv(cov)
    #ff.write('\n# inverse of covariance matrix \n')
    #for i in range(num):
    #    for j in range(num):
    #        ff.write(' %28.8f'%invC[i,j])
    #    ff.write('\n')

    L = np.linalg.cholesky(cov)
    B = np.linalg.solve(L, np.eye(num))
    invM = np.linalg.solve(np.transpose(L), B)
    ff.write('\n# inverse of covariance matrix \n')
    for i in range(num):
        for j in range(num):
            ff.write(' %28.8f'%invM[i,j])
        ff.write('\n')

    ff.close()
    print(corrs)
    #np.savetxt(roots+'_'+pinfo+'.txt', corrs, fmt='%16.8f')
    return

#===================================
if __name__ == "__main__":
    arg = []
    for ss in sys.argv[1:]:
        arg.append('%s'%ss)

    roots = arg[0]
    samples = loadMCSamples(roots, dist_settings={'ignore_rows':0.3})

    pname = samples.getParams()
    #la1 = np.pi*pname.DAstar*1000.0/pname.rstar
    #samples.addDerived(la1, name='lA1', label='l_A')
    la2 = np.pi/(pname.thetastar/100.0)
    samples.addDerived(la2, name='lA', label='l_A')
    R = np.sqrt(pname.omegam*pname.H0**2) *pname.DAstar*1000.0/(c_light/1000)
    samples.addDerived(R, name='R', label='R')

    samples.contours = np.array([0.68, 0.95, 0.99])
    samples.updateBaseStatistics()
    #triangleplot(samples)
    #outputcorrelation(roots,samples)
    outputdats(samples)
    
