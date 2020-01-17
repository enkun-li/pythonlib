#!/home/ekli/.local/anaconda3/bin/python
#-*- coding: utf-8 -*-  
#==================================
# File Name: my_gapp.py
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-06-11 11:09:46
#==================================

import ctypes
from ctypes import c_int, c_double, c_char, POINTER, c_char_p
import numpy as np
import matplotlib.pyplot as plt
import os

# Try to locate the .so file in the same directory as this file
_file = './libgapp.so'
_path = os.path.join(*(os.path.split(__file__)[:-1]+(_file,)))
_mod = ctypes.cdll.LoadLibrary(_path)

# int initial_gapp(N,filename,sigf,lenf)
initial_gapp = _mod.initial_gapp
initial_gapp.argtypes = (c_int, c_char_p)
initial_gapp.restype = c_int

# int initial_gapp_cov(char * filename);
initial_gapp_cov = _mod.initial_gapp_cov
initial_gapp_cov.argtypes = (c_char_p,)
initial_gapp_cov.restype = c_int

# int free_gapp()
free_gapp = _mod.free_gapp
free_gapp.restype = c_int

# int setup_gapp(double sigf, double lenf)
setup_gapp = _mod.setup_gapp
setup_gapp.argtypes(c_double, c_double)
setup_gapp.restype = c_int

# double rec_mu(x)
rec_mu = _mod.rec_mu
rec_mu.argtypes = (c_double,)
rec_mu.restype = c_double

# double rec_covariance(x,y)
rec_covariance = _mod.rec_covariance
rec_covariance.argtypes = (c_double, c_double)
rec_covariance.restype = c_double

# double rec_mu_over_mu0(x)
rec_mu_over_mu0 = _mod.rec_mu_over_mu0
rec_mu_over_mu0.argtypes = (c_double,)
rec_mu_over_mu0.restype = c_double

# double rec_cov_mu_over_mu0(x,y)
rec_cov_mu_over_mu0 = _mod.rec_cov_mu_over_mu0
rec_cov_mu_over_mu0.argtypes = (c_double, c_double)
rec_cov_mu_over_mu0.restype = c_double

# double integrate_one_over_mu(x)
integrate_one_over_mu = _mod.integrate_one_over_mu
integrate_one_over_mu.argtypes = (c_double,)
integrate_one_over_mu.restype = c_double

# double integrate_cov_one_over_mu(x,y)
integrate_cov_one_over_mu = _mod.integrate_cov_one_over_mu
integrate_cov_one_over_mu.argtypes = (c_double, c_double)
integrate_cov_one_over_mu.restype = c_double

# double rec_mu_int_one_over_mu(x)
rec_mu_int_one_over_mu = _mod.rec_mu_int_one_over_mu
rec_mu_int_one_over_mu.argtypes = (c_double,)
rec_mu_int_one_over_mu.restype = c_double

# double rec_cov_mu_int_one_over_mu(x,y)
rec_cov_mu_int_one_over_mu = _mod.rec_cov_mu_int_one_over_mu
rec_cov_mu_int_one_over_mu.argtypes = (c_double,c_double)
rec_cov_mu_int_one_over_mu.restype = c_double

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#double rec_fs8_Ii_Hz(double z, double gam);
rec_fs8_Ii_Hz = _mod.rec_fs8_Ii_Hz
rec_fs8_Ii_Hz.argtypes = (c_double, c_double)
rec_fs8_Ii_Hz.restype = c_double

#double rec_fs8_Mij_Hz(double x, double y, double gam);
rec_fs8_Mij_Hz = _mod.rec_fs8_Mij_Hz
rec_fs8_Mij_Hz.argtypes = (c_double, c_double, c_double)
rec_fs8_Mij_Hz.restype = c_double

#double rec_fs8_Nij_Hz(double zi, double zj, double gam);
rec_fs8_Nij_Hz = _mod.rec_fs8_Nij_Hz
rec_fs8_Nij_Hz.argtypes = (c_double, c_double, c_double)
rec_fs8_Nij_Hz.restype = c_double


#double rec_fs8_Ii_Ez(double z, double gam);
rec_fs8_Ii_Ez = _mod.rec_fs8_Ii_Ez
rec_fs8_Ii_Ez.argtypes = (c_double, c_double)
rec_fs8_Ii_Ez.restype = c_double

#double rec_fs8_Mij_Ez(double x, double y, double gam);
rec_fs8_Mij_Ez = _mod.rec_fs8_Mij_Ez
rec_fs8_Mij_Ez.argtypes = (c_double, c_double, c_double)
rec_fs8_Mij_Ez.restype = c_double

#double rec_fs8_Nij_Ez(double zi, double zj, double gam);
rec_fs8_Nij_Ez = _mod.rec_fs8_Nij_Ez
rec_fs8_Nij_Ez.argtypes = (c_double, c_double, c_double)
rec_fs8_Nij_Ez.restype = c_double


# double rec_fs8_with_Hz(x,sig8,ommh2,gam)
rec_fs8_with_Hz = _mod.rec_fs8_with_Hz
rec_fs8_with_Hz.argtypes = (c_double, c_double, c_double, c_double)
rec_fs8_with_Hz.restype = c_double

# double rec_cov_fs8_with_Hz(x,y,sig8,ommh2,gam)
rec_cov_fs8_with_Hz = _mod.rec_cov_fs8_with_Hz
rec_cov_fs8_with_Hz.argtypes = (c_double, c_double, c_double, c_double, c_double)
rec_cov_fs8_with_Hz.restype = c_double

# double rec_fs8_with_Hz(x,sig8,omm,gam)
rec_fs8_with_Ez = _mod.rec_fs8_with_Ez
rec_fs8_with_Ez.argtypes = (c_double, c_double, c_double, c_double)
rec_fs8_with_Ez.restype = c_double

# double rec_cov_fs8_with_Hz(x,y,sig8,omm,gam)
rec_cov_fs8_with_Ez = _mod.rec_cov_fs8_with_Ez
rec_cov_fs8_with_Ez.argtypes = (c_double, c_double, c_double, c_double, c_double)
rec_cov_fs8_with_Ez.restype = c_double

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# int initial_fs8_obs(int N, char *file_fz, char *file_cov);
initial_fs8_obs = _mod.initial_fs8_obs
initial_fs8_obs.argtypes = (c_int, c_char_p, c_char_p)
initial_fs8_obs.restype = c_int

# int free_fs8obs();
free_fs8obs = _mod.free_fs8obs
free_fs8obs.restype = c_int
 
# double fs8_loglikelihood_Hz(double sig8, double ommh2, double gam);
fs8_loglikelihood_Hz = _mod.fs8_loglikelihood_Hz
fs8_loglikelihood_Hz.argtypes = (c_double, c_double, c_double)
fs8_loglikelihood_Hz.restype = c_double

# double fs8_loglikelihood_Ez(double sig8, double omm, double gam);
fs8_loglikelihood_Ez = _mod.fs8_loglikelihood_Ez
fs8_loglikelihood_Ez.argtypes = (c_double, c_double, c_double)
fs8_loglikelihood_Ez.restype = c_double
 
# double fs8_loglikelihood_Hz_t(double sig8, double ommh2, double gam);
fs8_loglikelihood_Hz_t = _mod.fs8_loglikelihood_Hz_t
fs8_loglikelihood_Hz_t.argtypes = (c_double, c_double, c_double)
fs8_loglikelihood_Hz_t.restype = c_double

# double fs8_loglikelihood_Hz_not(double sig8, double ommh2, double gam);
fs8_loglikelihood_Hz_not = _mod.fs8_loglikelihood_Hz_not
fs8_loglikelihood_Hz_not.argtypes = (c_double, c_double, c_double)
fs8_loglikelihood_Hz_not.restype = c_double

# double fs8_loglikelihood_Ez_not(double sig8, double omm, double gam);
fs8_loglikelihood_Ez_not = _mod.fs8_loglikelihood_Ez_not
fs8_loglikelihood_Ez_not.argtypes = (c_double, c_double, c_double)
fs8_loglikelihood_Ez_not.restype = c_double


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# int initial_lcdm(omm,h0)
initial_lcdm = _mod.initial_lcdm
initial_lcdm.argtypes = (c_double, c_double)
initial_lcdm.restype = c_int

# double lcdm_hz(z)
lcdm_hz = _mod.lcdm_hz
lcdm_hz.argtypes = (c_double,)
lcdm_hz.restype = c_double

# double lcdm_comoving_distance(z)
lcdm_comoving_distance = _mod.lcdm_comoving_distance
lcdm_comoving_distance.argtypes = (c_double,)
lcdm_comoving_distance.restype = c_double

# double lcdm_growth(z,sig8,gam)
lcdm_growth = _mod.lcdm_growth
lcdm_growth.argtypes = (c_double, c_double, c_double)
lcdm_growth.restype = c_double




#===================================
if __name__ == "__main__":
    N = 51
    filename = bytes("/home/ekli/myworks/cosmodata/OHD_51.txt", "utf8")
    sigf=170.75260325
    lenf=2.61947963

    status = initial_gapp(N,filename,sigf,lenf)
    print(status)

    status = setup_gapp()
    print(status)
    
    h0 = rec_mu(0.0)
    omm = 0.3
    status = initial_lcdm(omm,h0)

    zz = np.linspace(0,2,41)

    for z in zz:
        print("%5.2f "%z, end=' ')
        print("%8.3f "%rec_mu(z), end=' ')
        cov = rec_covariance(z, z)
        print("%9.6f "%np.sqrt(cov), end=' ')
        print("%8.3f "%lcdm_hz(z), end=' ')
        
        print("")

    hz = np.array([rec_mu(z) for z in zz])
    sig = np.array([np.sqrt(rec_covariance(z,z)) for z in zz])

    plt.plot(zz, hz, '--b')
    plt.fill_between(zz, hz+sig, hz-sig, color='blue', alpha=0.5)

    lhz = np.array([lcdm_hz(z) for z in zz])
    plt.plot(zz, lhz, ':k')


    plt.show()

