#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: read_PLA_distancprior.py
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-04-04 12:14:44
#==================================
'''
    This file is used to read in datas from a txt file, in which some ignored
    columns will not read in.
    All the data will be read in to a dict.
'''


import numpy as np
import sys

def outputdat(root):
    ff = open(root, 'r')
    ss = ff.readlines()
    pars = {}
    pars['name'] = ss[0].split(sep='#')[1][:-2]
    nlines = len(ss)
    for i in range(nlines):
        if(ss[i][:7] == 'num_dat'):
            pars['num_dat'] = int( ss[i].split(sep='=')[-1] )
        elif(ss[i][:4] == '#Par'):
            par_line = i
        elif(ss[i][:4] == '#Cov'):
            cov_line = i
        elif(ss[i][:4] == '#Cor'):
            cor_line = i
        elif(ss[i][:4] == '#Inv'):
            inv_line = i
    ndat = pars['num_dat']
    cov = np.zeros((ndat, ndat))
    cor = np.zeros((ndat, ndat))
    invC = np.zeros((ndat, ndat))
    for i in range(ndat):
        par_dat = ss[par_line+i+1].split()
        pars[par_dat[0]] = [float(par_dat[1]), float(par_dat[2]),
                float(par_dat[3]), par_dat[4]]
        cov_dat = ss[cov_line+i+1].split()
        cov[i,:] = [float(cc) for cc in cov_dat]
        cor_dat = ss[cor_line+i+1].split()
        cor[i,:] = [float(cc) for cc in cor_dat]
        inv_dat = ss[inv_line+i+1].split()
        invC[i,:] = [float(cc) for cc in inv_dat]

    pars['cov'] = cov
    pars['corr'] = cor
    pars['invcov'] = invC
    return pars


#================================
if __name__ == '__main__':
    arg = sys.argv[1]
    pars = outputdat(arg)
    print(pars)
