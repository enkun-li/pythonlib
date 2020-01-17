#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: equations.pyx
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-04-18 10:20:57
#==================================

import sys
if(not sys.path[0][-10:-1] == 'pythonlib'):
    sys.path.insert(0, '/home/ekli/myworks/pythonlib/')
from constants import constants as const
import numpy as np
import yaml

cdef class LambdaGeneral:
    '''
    LambdaGeneral(params_info=None,feedback)
    =========================================================
    param params_info: a dict contains model parameters
    param feedback: 0 -> don't print anything
                    1 -> print some infomathons of model
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    example of params_file/info:
    '''
    cdef dict info
    cdef double omegab, omegac, omegag # baryon, CDM, photon
    cdef double omegav # DE

    def __cinit__(self,params_info=None,feedback=0):
        # read in parameters info
        if(params_info is not None):
            self.info = params_info
            if(feedback == 1):
                print("You have provided an info dict.")
        else:
            print("Error: you haven't provided parameters info. [LambdaGeneral]")
            sys.exit(0)

        if(feedback == 1):
            print()
        return

    cpdef double Omega_x(self, double a):
        '''
        Omega_x(a)
        =========================
        mod = 3: --> EDE model
        '''
        cdef ome, omx0, omx, omm, w0, omegax
        cdef int mod

        mod = self.info['mod']

        if(mod == 3):
            omm = self.info['omegam']
            omx = self.info['omegav']
            w0 = self.info['nps'][0]
            ome = self.info['nps'][1]
            omx0 = omm*omx/(1-omx)            
            omegax = (omx0 - ome*(1-a**(-3*w0)))/(omx0 +omm*a**(3*w0)) \
                    +ome*(1-a**(-3*w0))
        else:
            print('Error: no such model. [Omega_x]')
            omegax = 0.0
            sys.exit(0)
        return omegax


cdef class MassiveNu:
    '''
    MassiveNu(params_info=None,feedback)
    ========================================================
    Initial massive neutrinos
    '''
    cdef dict info

    def __cinit__(self,params_info=None,feedback=0):
        # read in parameters info
        if(params_info is not None):
            self.info = params_info
            if(feedback == 1):
                print("You have provided an info dict.")
        else:
            print("Error: you haven't provided parameters info. [LambdaGeneral]")
            sys.exit(0)

        return

    cpdef Nu_init(self, qmax):
        return

