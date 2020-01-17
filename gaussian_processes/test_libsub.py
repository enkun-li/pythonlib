#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: test_libsub.py
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-05-19 10:30:26
#==================================

from subroutines import *
import numpy as np

def numpy_kernel(x,y,sigf,lenf):
    return sigf**2 *np.exp( -(x-y)**2/2/lenf**2 )

for i in range(20):
    x,y = np.random.uniform(-1,1,2)

    k = kernel(x,y, 2, 1)
    r = numpy_kernel(x,y,2,1)
    print("%12.6f %12.6f %12.6f %12.6f"%(x,y,k,r))
