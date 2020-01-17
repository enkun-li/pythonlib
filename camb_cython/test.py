#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: test.py
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-04-15 19:47:49
#==================================

import numpy as np
import matplotlib.pyplot as plt
from modules import ModelParams
import yaml

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
infofile = './parameters.yaml'
CP = ModelParams(params_file=infofile, feedback=1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
infofile = './parameters.yaml'
#CP = ModelParams(params_file=infofile, feedback=1)
ff = open(infofile)
info = yaml.load(ff)
ff.close()

info['dark_energy_model'] = 3
CP = ModelParams(params_info=info, feedback=1)

zz = np.linspace(0,2,101)
hz = np.array([CP.Hofz(z) for z in zz])

plt.plot(zz, hz)

plt.show()
