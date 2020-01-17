#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: subroutines.py
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-05-19 10:27:43
#==================================

import os
import ctypes
from ctypes import c_int, c_double

# Try to locate the .so file in the same directory as this file
_file = './libsub.so'
_path = os.path.join(*(os.path.split(__file__)[:-1] + (_file,)))
_mod = ctypes.cdll.LoadLibrary(_path)

# double kernel(double, double)
kernel = _mod.kernel
kernel.argtypes = (c_double, c_double, c_double, c_double)
kernel.restype = c_double
