#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: pso.py
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-05-16 12:35:22
#==================================

import numpy as np
import matplotlib.pyplot as plt

class PSO:

    def __init__(self, nVar=5):
        ## Problem Definition
        self.nVar = nVar # Number of Unknown (Decision) Variables
        self.VarSize = n # Matrix Size of Decision Variables
        self.VarMin = # Lower Bound of Decision Variables
        self.VarMax = # Upper Bound of Decision Variables
        
        ## Parameters of PSO
        self.MaxIt = 100 # Maximum Number of Iterations
        self.nPop = 50 # Population Size (Swarm Size)
        self.w  = 1 # Intertia Coefficient
        self.c1 = 2 # Persinal Acceleration Coefficient
        self.c2 = 2 # Social Acceleration Coefficient
        
        ## Initialization
        empty_particle.Position = []
        empty_particle.Velocity = []
        empty_particle.Cost = []
        empty_particle.Best.Position = []
        empty_particle.Best.Cost = []

        par

        return

    ## Problem Definition
    def CostFunction(self, x, func):
        '''
        Cost Function
        '''
        return func(x)
    


## Initialization

## Main Loop of PSO

## Results


