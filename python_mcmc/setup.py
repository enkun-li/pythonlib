#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: setup.py
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2018-12-25 21:25:47
#==================================

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("mcmc", ["mcmc.pyx"]) ]

setup(
        name="mcmc pyx",
        cmdclass = {'build_ext': build_ext},
        ext_modules = ext_modules
)
