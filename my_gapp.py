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

# int initial_gapp_cov(int N, char * filename);
initial_gapp_cov = _mod.initial_gapp_cov
initial_gapp_cov.argtypes = (c_int, c_char_p)
initial_gapp_cov.restype = c_int

# int free_gapp()
free_gapp = _mod.free_gapp
free_gapp.restype = c_int

# int setup_gapp(double sigf, double lenf)
setup_gapp = _mod.setup_gapp
setup_gapp.argtypes = (c_double, c_double)
setup_gapp.restype = c_int

# double rec_mu(x)
rec_mu = _mod.rec_mu
rec_mu.argtypes = (c_double,)
rec_mu.restype = c_double

# double rec_covariance(x,y)
rec_covariance = _mod.rec_covariance
rec_covariance.argtypes = (c_double, c_double)
rec_covariance.restype = c_double

#double rec_distance_noM(double z);
rec_distance_noM = _mod.rec_distance_noM
rec_distance_noM.argtypes = (c_double,)
rec_distance_noM.restype = c_double

#double rec_cov_distance_noM(double zi, double zj);
rec_cov_distance_noM = _mod.rec_cov_distance_noM
rec_cov_distance_noM.argtypes = (c_double,c_double)
rec_cov_distance_noM.restype = c_double

#=========================================
#double rec_mu0_over_mu(double x, double * args);
#double rec_covxy_over_muxy(double x, double y);
#double rec_cov_mu0_over_mu(double x, double y);
#double rec_covmu0mu_over_fmu0mu(double x, double y, double * args);
#
#double int_mu0_over_mu(double x);
int_mu0_over_mu = _mod.int_mu0_over_mu
int_mu0_over_mu.argtypes = (c_double,)
int_mu0_over_mu.restype = c_double

#double cov_int_mu0_over_mu(double x, double y);
cov_int_mu0_over_mu = _mod.cov_int_mu0_over_mu
cov_int_mu0_over_mu.argtypes = (c_double, c_double)
cov_int_mu0_over_mu.restype = c_double

#int return_Pantheon_dz(void);

#=========================================
# double rec_mu_over_mu0(x)
##rec_mu_over_mu0 = _mod.rec_mu_over_mu0
##rec_mu_over_mu0.argtypes = (c_double,)
##rec_mu_over_mu0.restype = c_double
##
### double rec_cov_mu_over_mu0(x,y)
##rec_cov_mu_over_mu0 = _mod.rec_cov_mu_over_mu0
##rec_cov_mu_over_mu0.argtypes = (c_double, c_double)
##rec_cov_mu_over_mu0.restype = c_double
##
### double integrate_one_over_mu(x)
##integrate_one_over_mu = _mod.integrate_one_over_mu
##integrate_one_over_mu.argtypes = (c_double,)
##integrate_one_over_mu.restype = c_double
##
### double integrate_cov_one_over_mu(x,y)
##integrate_cov_one_over_mu = _mod.integrate_cov_one_over_mu
##integrate_cov_one_over_mu.argtypes = (c_double, c_double)
##integrate_cov_one_over_mu.restype = c_double
##
### double rec_mu_int_one_over_mu(x)
##rec_mu_int_one_over_mu = _mod.rec_mu_int_one_over_mu
##rec_mu_int_one_over_mu.argtypes = (c_double,)
##rec_mu_int_one_over_mu.restype = c_double
##
### double rec_cov_mu_int_one_over_mu(x,y)
##rec_cov_mu_int_one_over_mu = _mod.rec_cov_mu_int_one_over_mu
##rec_cov_mu_int_one_over_mu.argtypes = (c_double,c_double)
##rec_cov_mu_int_one_over_mu.restype = c_double
##
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###double rec_fs8_Ii_Hz(double z, double gam);
##rec_fs8_Ii_Hz = _mod.rec_fs8_Ii_Hz
##rec_fs8_Ii_Hz.argtypes = (c_double, c_double)
##rec_fs8_Ii_Hz.restype = c_double
##
###double rec_fs8_Mij_Hz(double x, double y, double gam);
##rec_fs8_Mij_Hz = _mod.rec_fs8_Mij_Hz
##rec_fs8_Mij_Hz.argtypes = (c_double, c_double, c_double)
##rec_fs8_Mij_Hz.restype = c_double
##
###double rec_fs8_Nij_Hz(double zi, double zj, double gam);
##rec_fs8_Nij_Hz = _mod.rec_fs8_Nij_Hz
##rec_fs8_Nij_Hz.argtypes = (c_double, c_double, c_double)
##rec_fs8_Nij_Hz.restype = c_double
##
##
###double rec_fs8_Ii_Ez(double z, double gam);
##rec_fs8_Ii_Ez = _mod.rec_fs8_Ii_Ez
##rec_fs8_Ii_Ez.argtypes = (c_double, c_double)
##rec_fs8_Ii_Ez.restype = c_double
##
###double rec_fs8_Mij_Ez(double x, double y, double gam);
##rec_fs8_Mij_Ez = _mod.rec_fs8_Mij_Ez
##rec_fs8_Mij_Ez.argtypes = (c_double, c_double, c_double)
##rec_fs8_Mij_Ez.restype = c_double
##
###double rec_fs8_Nij_Ez(double zi, double zj, double gam);
##rec_fs8_Nij_Ez = _mod.rec_fs8_Nij_Ez
##rec_fs8_Nij_Ez.argtypes = (c_double, c_double, c_double)
##rec_fs8_Nij_Ez.restype = c_double
##
##
### double rec_fs8_with_Hz(x,sig8,ommh2,gam)
##rec_fs8_with_Hz = _mod.rec_fs8_with_Hz
##rec_fs8_with_Hz.argtypes = (c_double, c_double, c_double, c_double)
##rec_fs8_with_Hz.restype = c_double
##
### double rec_cov_fs8_with_Hz(x,y,sig8,ommh2,gam)
##rec_cov_fs8_with_Hz = _mod.rec_cov_fs8_with_Hz
##rec_cov_fs8_with_Hz.argtypes = (c_double, c_double, c_double, c_double, c_double)
##rec_cov_fs8_with_Hz.restype = c_double
##
### double rec_fs8_with_Hz(x,sig8,omm,gam)
##rec_fs8_with_Ez = _mod.rec_fs8_with_Ez
##rec_fs8_with_Ez.argtypes = (c_double, c_double, c_double, c_double)
##rec_fs8_with_Ez.restype = c_double
##
### double rec_cov_fs8_with_Hz(x,y,sig8,omm,gam)
##rec_cov_fs8_with_Ez = _mod.rec_cov_fs8_with_Ez
##rec_cov_fs8_with_Ez.argtypes = (c_double, c_double, c_double, c_double, c_double)
##rec_cov_fs8_with_Ez.restype = c_double
##
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### int initial_fs8_obs(int N, char *file_fz, char *file_cov);
##initial_fs8_obs = _mod.initial_fs8_obs
##initial_fs8_obs.argtypes = (c_int, c_char_p, c_char_p)
##initial_fs8_obs.restype = c_int
##
### int free_fs8obs();
##free_fs8obs = _mod.free_fs8obs
##free_fs8obs.restype = c_int
## 
### double fs8_loglikelihood_Hz(double sig8, double ommh2, double gam);
##fs8_loglikelihood_Hz = _mod.fs8_loglikelihood_Hz
##fs8_loglikelihood_Hz.argtypes = (c_double, c_double, c_double)
##fs8_loglikelihood_Hz.restype = c_double
##
### double fs8_loglikelihood_Ez(double sig8, double omm, double gam);
##fs8_loglikelihood_Ez = _mod.fs8_loglikelihood_Ez
##fs8_loglikelihood_Ez.argtypes = (c_double, c_double, c_double)
##fs8_loglikelihood_Ez.restype = c_double
## 
### double fs8_loglikelihood_Hz_t(double sig8, double ommh2, double gam);
##fs8_loglikelihood_Hz_t = _mod.fs8_loglikelihood_Hz_t
##fs8_loglikelihood_Hz_t.argtypes = (c_double, c_double, c_double)
##fs8_loglikelihood_Hz_t.restype = c_double
##
### double fs8_loglikelihood_Hz_not(double sig8, double ommh2, double gam);
##fs8_loglikelihood_Hz_not = _mod.fs8_loglikelihood_Hz_not
##fs8_loglikelihood_Hz_not.argtypes = (c_double, c_double, c_double)
##fs8_loglikelihood_Hz_not.restype = c_double
##
### double fs8_loglikelihood_Ez_not(double sig8, double omm, double gam);
##fs8_loglikelihood_Ez_not = _mod.fs8_loglikelihood_Ez_not
##fs8_loglikelihood_Ez_not.argtypes = (c_double, c_double, c_double)
##fs8_loglikelihood_Ez_not.restype = c_double
##
##
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# int initial_lcdm(omm,omk,h0)
initial_lcdm = _mod.initial_lcdm
initial_lcdm.argtypes = (c_double, c_double, c_double)
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

#double lcdm_como_ang_diam_distance(double z); // in Mpc
lcdm_como_ang_diam_distance = _mod.lcdm_como_ang_diam_distance
lcdm_como_ang_diam_distance.argtypes = (c_double,)
lcdm_como_ang_diam_distance.restype = c_double

#double lcdm_luminosity_distance(double z);
lcdm_luminosity_distance = _mod.lcdm_luminosity_distance
lcdm_luminosity_distance.argtypes = (c_double,)
lcdm_luminosity_distance.restype = c_double

#double lcdm_angulardiameter_distance(double z);
lcdm_angulardiameter_distance = _mod.lcdm_angulardiameter_distance
lcdm_angulardiameter_distance.argtypes = (c_double,)
lcdm_angulardiameter_distance.restype = c_double

#double lcdm_mu_SN(double z);
lcdm_mu_SN = _mod.lcdm_mu_SN
lcdm_mu_SN.argtypes = (c_double,)
lcdm_mu_SN.restype = c_double

#====================================================
#int initial_Pantheon(int N, char * file_dat, char * file_cov);
initial_Pantheon = _mod.initial_Pantheon
initial_Pantheon.argtypes = (c_int, c_char_p, c_char_p)
initial_Pantheon.restype = c_int

#int free_Pantheon(void);
free_Pantheon = _mod.free_Pantheon
free_Pantheon.restype = c_int

#double sn_loglikelihood(double (*func)(double));
#sn_loglikelihood
#double lcdm_pantheon(double h0, double omm, double omk);
lcdm_pantheon = _mod.lcdm_pantheon
lcdm_pantheon.argtypes = (c_double, c_double, c_double)
lcdm_pantheon.restype = c_double

#==================================================
#int initial_SL_data(int N, char *filename);
#double sl_loglikelihood(double (*func)(double, double), double fE);
#int free_stronglens(void);

#double lcdm_stronglens(double h0, double omm, double omk, double fE);
lcdm_stronglens = _mod.lcdm_stronglens
lcdm_stronglens.argtypes = (c_double, c_double, c_double, c_double)
lcdm_stronglens.restype = c_double

#double cosmography_stronglens(double omk, double a1, double a2,  double fE);
cosmography_stronglens = _mod.cosmography_stronglens
cosmography_stronglens.argtypes = (c_double, c_double, c_double, c_double)
cosmography_stronglens.restype = c_double

#=======================================
#int initial_SL_loglike(int N, char * file_cov);
initial_SL_loglike = _mod.initial_SL_loglike
initial_SL_loglike.argtypes = (c_int, c_char_p)
initial_SL_loglike.restype = c_int

#int free_SL(void);
free_SL = _mod.free_SL
free_SL.restype = c_int

#double SL_loglike(double MB, double omk, double fE);
SL_loglike = _mod.SL_loglike
SL_loglike.argtypes = (c_double, c_double, c_double)
SL_loglike.restype = c_double

#int initial_TDSL_loglike(int N, char * file_dat, char * file_cov);
initial_TDSL_loglike = _mod.initial_TDSL_loglike
initial_TDSL_loglike.argtypes = (c_int, c_char_p, c_char_p)
initial_TDSL_loglike.restype = c_int

#int free_TDSL(void);
free_TDSL = _mod.free_SL
free_TDSL.restype = c_int

#double TDSL_loglike(double MB, double H0, double omk);
TDSL_loglike = _mod.TDSL_loglike
TDSL_loglike.argtypes = (c_double, c_double, c_double)
TDSL_loglike.restype = c_double

#double return_SL_loglike(double MB, double omk, double fE);
return_SL_loglike = _mod.return_SL_loglike
return_SL_loglike.argtypes = (c_double, c_double, c_double)
return_SL_loglike.restype = c_double

#double return_TDSL_loglike(double MB, double omk, double H0);
return_TDSL_loglike = _mod.return_TDSL_loglike
return_TDSL_loglike.argtypes = (c_double, c_double, c_double)
return_TDSL_loglike.restype = c_double

#double margin_of_MB_SL_loglike(double omk, double fE);
margin_of_MB_SL_loglike = _mod.margin_of_MB_SL_loglike
margin_of_MB_SL_loglike.argtypes = (c_double, c_double)
margin_of_MB_SL_loglike.restype = c_double

#double margin_of_MB_TDSL_loglike(double omk, double H0);
margin_of_MB_TDSL_loglike = _mod.margin_of_MB_TDSL_loglike
margin_of_MB_TDSL_loglike.argtypes = (c_double, c_double)
margin_of_MB_TDSL_loglike.restype = c_double

#===========================================================
##ifndef __COSMO_GP__
##define __COSMO_GP__
#
#int initial_cosmo_gp(int num);
initial_cosmo_gp = _mod.initial_cosmo_gp
initial_cosmo_gp.argtypes = (c_int,)
initial_cosmo_gp.restype = c_int

#int setup_cosmo_gp(double omegak, int ohdtype);
setup_cosmo_gp = _mod.setup_cosmo_gp
setup_cosmo_gp.argtypes = (c_double, c_double)
setup_cosmo_gp.restype = c_int

#double gp_hz(double z);
gp_hz = _mod.gp_hz
gp_hz.argtypes = (c_double,)
gp_hz.restype = c_double

#double gp_cov_hz(double zi, double zj);
gp_cov_hz = _mod.gp_cov_hz
gp_cov_hz.argtypes = (c_double,c_double)
gp_cov_hz.restype = c_double

#double gp_dimles_comdistance(double z);
gp_dimles_comdistance = _mod.gp_dimles_comdistance
gp_dimles_comdistance.argtypes = (c_double,)
gp_dimles_comdistance.restype = c_double

#double gp_cov_dimles_comdistance(double zi, double zj);
gp_cov_dimles_comdistance = _mod.gp_cov_dimles_comdistance
gp_cov_dimles_comdistance.argtypes = (c_double,c_double)
gp_cov_dimles_comdistance.restype = c_double

#double gp_dimles_comangdistance(double z);
#double gp_cov_dimles_comangdistance(double zi, double zj);
#
#double gp_com_distance(double z);
#double gp_ang_distance(double z);
#double gp_lum_distance(double z);
#
#double gp_DistanceModulus_star(double z);
#double gp_cov_DistanceModulus_star(double zi, double zj);
#
#double gp_DistanceModulus(double z);
#
#double gp_DistanceSumRole(double zl, double zs);
#double gp_cov_DistanceSumRole(double zli, double zsi, double zlj, double zsj);
#
#int GP_OHD_Pantheon_output_dc(int num);
#
#int GP_OHD_initial_Pantheon(void);
#double GP_OHD_sn_loglikelihood(void);
#double GP_OHD_rofchi(double chi, double omk);
#double return_gp_ohd_loglike(double omk);
return_gp_ohd_loglike = _mod.return_gp_ohd_loglike
return_gp_ohd_loglike.argtypes = (c_double,)
return_gp_ohd_loglike.restype = c_double

#endif /* __COSMO_GP__ */
###===================================
##if __name__ == "__main__":
##    N = 51
##    filename = bytes("/home/ekli/myworks/cosmodata/OHD_51.txt", "utf8")
##    sigf=170.75260325
##    lenf=2.61947963
##
##    status = initial_gapp(N,filename,sigf,lenf)
##    print(status)
##
##    status = setup_gapp()
##    print(status)
##    
##    h0 = rec_mu(0.0)
##    omm = 0.3
##    status = initial_lcdm(omm,h0)
##
##    zz = np.linspace(0,2,41)
##
##    for z in zz:
##        print("%5.2f "%z, end=' ')
##        print("%8.3f "%rec_mu(z), end=' ')
##        cov = rec_covariance(z, z)
##        print("%9.6f "%np.sqrt(cov), end=' ')
##        print("%8.3f "%lcdm_hz(z), end=' ')
##        
##        print("")
##
##    hz = np.array([rec_mu(z) for z in zz])
##    sig = np.array([np.sqrt(rec_covariance(z,z)) for z in zz])
##
##    plt.plot(zz, hz, '--b')
##    plt.fill_between(zz, hz+sig, hz-sig, color='blue', alpha=0.5)
##
##    lhz = np.array([lcdm_hz(z) for z in zz])
##    plt.plot(zz, lhz, ':k')
##
##
##    plt.show()
##
