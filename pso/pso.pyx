#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: pso.pyx
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-03-22 10:01:29
#==================================

cimport cython
import numpy as np

#====================================================
# Particle Swarm Optimization
#
cdef class PSO:
    '''
    PSO()
    =====================================
    :: Particle Swarm Optimization
    param w: inertia weight
    param c1: social cognition learning factor
    param c2: self-knowledge learning factor
    param maxgen: number of evolution/iteration
    param sizepop: population size / number of particles
    '''
    cdef double w, c1, c2 # 
    cdef int maxgen, sizepop #
