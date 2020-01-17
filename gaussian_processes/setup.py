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

ext_modules = [
        Extension("gaussian_process", ["gaussian_process.pyx"]),
        Extension("gpsn_likelihood", ["gpsn_likelihood.pyx"]),
        Extension("gp_likelihood", 
            ["gp_likelihood.pyx"],
            extra_compile_args=['-fopenmp'],
            extra_link_args=['-fopenmp'],
            ),
        Extension("gp_c_extern", ["gp_c_extern.pyx"]),
        Extension("general_fsigma", ["general_fsigma.pyx"]),
        #Extension("general_fsigma_noH0", ["general_fsigma_noH0.pyx"]),
        ]

setup(
        name="gaussian_process",
        cmdclass = {'build_ext': build_ext},
        ext_modules = ext_modules
)
