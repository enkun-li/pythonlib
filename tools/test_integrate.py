#!/home/ekli/.local/anaconda3/bin/python
#-*- coding: utf-8 -*-  
#==================================
# File Name: test_integrate.py
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-06-10 09:53:49
#==================================

import sys
from integratefunc import *

def print_x_w(x,w,n):
    print("%5s %12s %12s"%('#', 'x[i]', 'w[i]'))
    for i in range(n+1):
        print('%5d %12.9f %12.9f'%(i, x[i], w[i]))
    return

def func(x):
    return x*np.exp(-x)

def actual_func(x1,x2):
    return (1+x1)*np.exp(-x1) -(1+x2)*np.exp(-x2)

def func_2d(x, y):
    return x*x + y*y +x*y

def actual_func_2d(a, b, c, d):
    return 1/12*(a-b)*(c-d)*(4*a*a+4*b*b+3*b*(c+d)+4*(c*c+c*d+d*d)+a*(4*b+3*(c+d)))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__=="__main__":
    arg = sys.argv[1:]
    x1 = float(arg[0])
    x2 = float(arg[1])
    n = int(arg[2])
    
    x, w = gauleg(x1,x2)
    print_x_w(x,w, n)

    ss = qgausleg(func, x1, x2)
    rs = actual_func(x1, x2)
    print("="*60)
    print("Actual value is %12.6f"%rs)
    print("Integrate value %12.6f"%ss)
    
    ss = qgausleg2d_fix(func_2d, 0.0, x1, 0.0, x2)
    rs = actual_func_2d(0.0, x1, 0.0, x2)
    
    print("="*60)
    print("Actual value is %12.6f"%rs)
    print("Integrate value %12.6f"%ss)
