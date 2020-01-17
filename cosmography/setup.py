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

ext_modules = [ \
        #Extension("darksector", ["darksector.pyx"]), \
        #Extension("cosmography", ["cosmography.pyx"]), \
        #Extension("likelihood", ["likelihood.pyx"]), \
        Extension("Model.lcdm_omk", ["Model_lcdm_omk.pyx"]),
        Extension("Model.lcdm_cgr", ["Model_lcdm_cgr.pyx"]),
        Extension("Model.cosm_tay", ["Model_cosm_tay.pyx"]),
        Extension("Model.cosm_pade", ["Model_cosm_pade.pyx"]),
        Extension("Likelihood.HST", ["Likelihood_HST.pyx"]),
        Extension("Likelihood.OHD", ["Likelihood_OHD.pyx"]),
        Extension("Likelihood.JLA", ["Likelihood_JLA.pyx"]),
        Extension("Likelihood.SDF", ["Likelihood_6dF.pyx"]),
        Extension("Likelihood.MGS", ["Likelihood_MGS.pyx"]),
        Extension("Likelihood.DR12Consensus", ["Likelihood_DR12Consensus.pyx"]),
        Extension("Likelihood.DR11quasars", ["Likelihood_DR11quasars.pyx"]),
        Extension("Likelihood.DR12Lya", ["Likelihood_DR12Lya.pyx"]),
        Extension("Likelihood.DR14quasars", ["Likelihood_DR14quasars.pyx"]),
        Extension("Likelihood.DR11Lya_xcro", ["Likelihood_DR11Lya_xcro.pyx"]),
        Extension("Likelihood.DR11Lya_auto", ["Likelihood_DR11Lya_auto.pyx"]),
        #Extension("neutrinos", ["neutrinos.pyx"])
        ]

setup(
        name="cosmography pyx",
        cmdclass = {'build_ext': build_ext},
        ext_modules = ext_modules
)
