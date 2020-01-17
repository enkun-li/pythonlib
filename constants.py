#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: constants.py
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-04-09 11:47:51
#==================================
'''
    This file is some general constants used in cosmology.
    It is a dict, you can use it just as: 
        `import constants as const`
    and use some constant in the form (such as `const_pi`)
        `const['const_pi']`
    enjoy!
'''

import numpy as np

constants = {}

constants['const_pi'] = 3.1415926535897932384626433832795
constants['const_twopi'] =2.*constants['const_pi']
constants['const_fourpi'] =4.*constants['const_pi']
constants['const_sqrt6'] =2.4494897427831780981972840747059

constants['c'] = 2.99792458e8
constants['h_P'] = 6.62606896e-34
constants['hbar'] = constants['h_P']/constants['const_twopi']

constants['G'] =6.6738e-11 #data book 2012, last digit +/-8
constants['sigma_thomson'] = 6.6524616e-29
constants['sigma_boltz'] = 5.6704e-8
constants['k_B'] = 1.3806504e-23
constants['eV'] = 1.60217646e-19


constants['m_p'] = 1.672621637e-27 # 1.672623e-27
constants['m_H'] = 1.673575e-27 #av. H atom
constants['m_e'] = 9.10938215e-31
constants['mass_ratio_He_H'] = 3.9715


constants['Gyr'] =3.1556926e16
constants['Mpc'] = 3.085678e22 #seem to be different definitions of this?
constants['Mpc_in_sec'] = constants['Mpc']/constants['c'] # Mpc/c = 1.029272d14 in SI units

constants['barssc0'] = constants['k_B'] / constants['m_p'] / constants['c']**2
constants['kappa'] =8.*constants['const_pi']*constants['G']
constants['a_rad'] = 8.*constants['const_pi']**5*constants['k_B']**4/15/constants['c']**3/constants['h_P']**3
#7.565914e-16 #radiation constant for u=aT^4


constants['Compton_CT'] = constants['Mpc_in_sec']*(8.0/3.0)*(constants['sigma_thomson']/(constants['m_e']*constants['c']))*constants['a_rad']
#Compton_CT is CT in Mpc units, (8./3.)*(sigma_T/(m_e*C))*a_R in Mpc
#Used to get evolution of matter temperature

#For 21cm
constants['f_21cm'] = 1420.40575e6
constants['l_21cm'] = constants['c']/constants['f_21cm']
constants['T_21cm'] = constants['h_P']*constants['f_21cm']/constants['k_B']
constants['A10'] = 2.869e-15
constants['B10'] = constants['l_21cm']**3/2/constants['h_P']/constants['c']*constants['A10']

constants['line21_const'] = 3*constants['l_21cm']**2*constants['c']*constants['h_P']/32/constants['const_pi']/constants['k_B']*constants['A10'] * constants['Mpc_in_sec'] * 1000
#1000 to get in MiliKelvin
constants['COBE_CMBTemp'] = 2.7255 #(Fixsen 2009) used as default value
constants['default_nnu'] = 3.046
#Neutrino mass splittings
constants['delta_mnu21'] = 7.54e-5 #eV^2 Particle Data Group 2015 (-0.22, + 0.26)
constants['delta_mnu31'] = 2.46e-3 #eV^2 Particle Data Group 2015 (+- 0.06)
#Round up to 0.06, so consistent with cosmomc's 1 neutrino default
constants['mnu_min_normal'] = 0.06 # sqrt(delta_mnu31)+sqrt(delta_mnu21)

constants['zeta3']  = 1.2020569031595942853997
constants['zeta5']  = 1.0369277551433699263313
constants['zeta7']  = 1.0083492773819228268397
