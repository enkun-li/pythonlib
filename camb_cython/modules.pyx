#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: modules.pyx
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-04-18 19:35:32
#==================================

import sys
if(not sys.path[0][-10:-1] == 'pythonlib'):
    sys.path.insert(0, '/home/ekli/myworks/pythonlib/')
from constants import constants as const
from equations import LambdaGeneral
import numpy as np
import yaml
from scipy import integrate

#===============================================
cdef class Initial_MP:
    '''
    LambdaGeneral(params_file=None,params_info=None,feedback)
    =========================================================
    param params_file: a yaml file contains model parameters
                       infomation
    param params_info: a dict contains model parameters, you 
                       can only use one params_file or 
                       params_info, if two are given, 
                       params_info will be the choosen one
    param feedback: 0 -> don't print anything
                    1 -> print some infomathons of model
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    example of params_file/info:
    '''
    cdef dict info
    cdef double omegab, omegac, omegag # baryon, CDM, photon
    cdef double omega_r # per neutrino species
    cdef double omegan # neutrinos
    cdef double omeganunomass # neutrinos no mass
    cdef double[:] omeganumass # neutrinos have mass
    cdef double omegav # DE
    cdef double omegak # spatial curvature
    cdef double omegam, omegar # matter (b,c,nu_massive) & radiation (g,nu_nomass)
    cdef double H0 # Hubble constant
    cdef double w_lam, wa_ppf # equation of state for DE
    cdef double Tcmb, Tnu # temperature of CMB, CNB (cosmic neutrino background)
    cdef double r, curv # for spatial curvature
    cdef int flat, Ksign # {closed,flat,open} -> {1,0,-1} -> {-1,0,1}
    cdef int mod # for different dark energy model
    cdef double[:] nps # array of new model parameters
    cdef double N_eff # effective number of neutrinos species
    cdef double[:] nu_mass, nu_mass_degeneracies # for massive neutrinos
    cdef int num_nu_mass # number of all massive neutrinos species
    cdef double num_nu_massless, nu_massless_degeneracy # for massless neutrinos


    def __cinit__(self,params_file=None,params_info=None,feedback=0):
        # read in parameters info
        if(params_file is not None):
            ff = open(params_file)
            self.info = yaml.load(ff)
            ff.close()
            if(feedback == 1):
                print("The parameters file is: %s"%params_file)
        elif(params_info is not None):
            self.info = params_info
            if(feedback == 1):
                print("You have provided an info dict.")
        else:
            print("Error: you haven't provided parameters info. [LambdaGeneral]")
            sys.exit(0)

        self.setup_parameters()
        if(feedback == 1):
            print("="*40)
            flat = ['closed', 'flat', 'open']
            print("%25s : %12s"%('The universe is', flat[self.flat+1]))
            demodel = self.info['DE_models']
            print("%25s : %12s"%('cosmological model', demodel[self.mod]))
            print("%25s = %12.6f"%('Om_b h^2', self.omegab*(self.H0/100)**2))
            print("%25s = %12.6f"%('Om_c h^2', self.omegac*(self.H0/100)**2))
            print("%25s = %12.6f"%('Om_nu h^2', self.omegan*(self.H0/100)**2))
            print("%25s = %12.6f"%('Om_x', self.omegav))
            print("%25s = %12.6f"%('Om_k', self.omegak))
            print("%25s = %12.6f"%('Om_g', self.omegag))
            print("~"*40)
            print("%25s = %12.6f"%('temperature of T_nu', self.Tnu) )
            print('%25s : %12d'%('Number of massless nu', self.num_nu_massless))
            print('%25s : %12d'%('degeneracy of massless nu', self.num_nu_massless))
            print('%25s : %12d'%('Number of massive nu', self.num_nu_mass))
            for i in range(self.num_nu_mass):
                print('%24s%d = %12.6f (%12.6f)'%('effective mass of massive nu_',\
                        i, self.nu_mass[i], self.nu_mass[i]*self.Tnu))
            print("%25s = %12.6f"%('total Om_nu,mass', \
                    np.sum(self.nu_mass)/94.15/(self.H0/100)**2*self.Tnu))
            print("%25s = %12.6f"%('Om_nu,nomass', self.omeganunomass))
            print("%25s = %12.6f"%('Om_r (per neutrino)', self.omega_r))
            print("%25s = %12.6f"%('Om_m (1-Om_k-Om_x)', 1-self.omegak-self.omegav))
            print("%25s = %12.6f"%('Om_m (Om_b+Om_c+Om_nu,massive)', self.omegab+self.omegac+self.omegan))
            print("%25s = %12s"%('Number of params user add', len(self.info['nps']) ))          
                
            return

    def print_input_info(self):
        '''
        print_input_info()
        ======================
        Output your input_info, just a test.
        '''
        print(self.info)
        return

    def setup_parameters(self):
        '''
        setup_parameters()
        ========================
        Setup the comman parameters the model need.
        '''
        cdef int num_nps
        cdef double fractional_numbel, neff_i
        cdef int i, actual_massless

        # read in base model parameters
        use_physical = self.info['use_physical']
        if(use_physical):
            self.H0 = self.info['h']*100
            self.omegab = self.info['ombh2']/(self.H0/100)**2
            self.omegac = self.info['omch2']/(self.H0/100)**2
            self.omegan = self.info['omnuh2']/(self.H0/100)**2
            self.omegak = self.info['omegak']
            self.omegav = 1-self.omegak-self.omegab-self.omegac-self.omegan
        else:
            self.H0 = self.info['H0']
            self.omegab = self.info['omegab']
            self.omegac = self.info['omegac']
            self.omegan = self.info['omegan']
            self.omegav = self.info['omegav']
            self.omegak = 1-self.omegab-self.omegac-self.omegav-self.omegan

        # temperature of CMB
        self.Tcmb= self.info['Tcmb']
        # energy density of photon \rho_g = \pi^2/15 * T_g^4
        # \rho_cr = 3c^2 H_0^2/(8\pi G) = 300 c^2/(8\pi G) h^2
        self.omegag = 2.38e-5*(self.Tcmb/2.7)**4/(self.H0/100.0)**2
        # energy density of single relativistic neutrino
        # self.omega_r = 7.0/8*(4.0/11)**(4.0/3)*self.omegag
        self.omega_r = 0.6193835936188336*self.omegag

        # effective of neutrinos species number
        self.N_eff = self.info['N_eff']
        # array of massive neutrinos mass
        self.nu_mass = np.array(self.info['nu_mass'])
        # number of massive neutrinos
        self.num_nu_mass = self.nu_mass.shape[0]
        # each mass corresponding to one degeneracy
        self.nu_mass_degeneracies = np.ones(self.num_nu_mass)
        # number of massless neutrino
        self.num_nu_massless = self.N_eff - np.float64(self.num_nu_mass)
        # number or degeneracy of massless neutrinos, 
        # all the same, so use one degeneracy
        self.nu_massless_degeneracy = self.num_nu_massless
        
        # whether share delta Neff: N_eff = num_nu_mass + num_nu_massless
        if(self.info['share_delta_neff']):
            fractional_numbel = self.N_eff
            actual_massless = int(self.num_nu_massless+1e-6)
            neff_i = np.float64(fractional_numbel)/ \
                    np.float64(actual_massless +self.num_nu_mass)
            self.num_nu_massless = neff_i*actual_massless
            self.nu_mass_degeneracies = np.array([self.nu_mass_degeneracies[i] \
                    *neff_i for i in range(self.num_nu_mass)])
            #print('share_delta_neff', fractional_numbel, actual_massless, neff_i)

        self.omeganunomass = self.omega_r*self.num_nu_massless
        self.omeganumass = np.zeros(self.num_nu_mass)
        # temperature of neutrinos
        self.Tnu = (4.0/11.0)**(1.0/3.0) *self.Tcmb
        for i in range(self.num_nu_mass):
            self.omeganumass[i]= self.omega_r * \
                    self.nu_mass_degeneracies[i]
            self.nu_mass[i] = self.nu_mass[i] /self.Tnu
        self.w_lam= self.info['w_lam']
        self.wa_ppf= self.info['wa_ppf']
        
        # flat, closed, open universe
        if(self.omegak > 1e-8):
            #open
            self.curv = self.H0**2 *self.omegak/(const['c']/1000)
            self.r = 1.0/np.sqrt(self.curv)
            self.Ksign = -1
            self.flat = 1
        elif(self.omegak < -1e-8):
            #closed
            self.curv = -self.H0**2 *self.omegak/(const['c']/1000)
            self.r = 1.0/np.sqrt(self.curv)
            self.Ksign = 1
            self.flat = -1
        #if(np.abs(self.omegak) < 1e-8 ):
        else:
            #flat
            self.curv = 0.0
            self.r = 1.0
            self.Ksign = 0
            self.flat = 0
            
        # determine the dark energy model
        self.mod = self.info['dark_energy_model']

        # the general matter energy density and radiation
        self.omegam = self.omegab +self.omegac
        self.omegar = self.omegag +self.omeganunomass

        # determine the new model parameters
        num_nps = len(self.info['nps'])
        self.nps = np.zeros(num_nps)
        for i in range(num_nps):
            self.nps[i] = self.info['nps'][i]

        return


#===============================================
cdef class ModelParams(Initial_MP):
    cdef dict dmde_info

    def __cinit__(self,params_file=None,params_info=None,feedback=0):
        cdef double zz, aa

        super(Initial_MP, self).__init__(params_file,params_info,feedback)

        self.dmde_info = {
                'omegam': self.omegam,
                'omegav': self.omegav,
                'nps': self.nps,
                'mod': self.mod}

        if(feedback == 1):
            print('~'*40)
            #print('ModelParams: Omega_m', self.omegab + self.omegac)
            zz = self.redshift_star()
            aa = 1/(1.0+zz)
            print('%25s = %12.6f'%('z_star', zz) )
            print('%25s = %12.6f'%('D_M(zstar) [Gpc]', self.AngularDiameterDistance(zz)*(1+zz)/1000) )
            print('%25s = %12.6f'%('r_s(zstar)', self.sound_horizon(1e-8, aa) ))
            print('%25s = %12.6f'%('100 theta_*', self.theta_star()*100) )
            print('~'*40)
            zz = self.redshift_drag()
            aa = 1/(1.0+zz)
            print('%25s = %12.6f'%('z_drag', zz) )
            print('%25s = %12.6f'%('D_M(zdrag) [Gpc]', self.AngularDiameterDistance(zz)*(1+zz)/1000) )
            print('%25s = %12.6f'%('r_s(zdrag)', self.sound_horizon(1e-8, aa) ))
            print('%25s = %12.6f'%('100 theta_d', self.theta_drag()*100) )
        return
    
    cpdef double Hofz(self, double z):
        '''
        Hofz(z)
        ===========================
        Hubble parameter in redshift z [H/c light speed in km/s]
        param z: redshift
        '''
        cdef double a, Ez2, hz
        cdef double omm, omx

        a = 1.0/(1+z)

        Ez2 = self.omegak/a**2 + self.omegab/a**3 +(self.omegag+self.omeganunomass)/a**4
        if(self.mod == 0 or self.mod == 1 or self.mod == 2):
            omm = self.omegac/a**3
            omx = self.omegav/a**(3*(1+self.w_lam+self.wa_ppf)) \
                    * np.exp(-3*self.wa_ppf*(1-a))
            Ez2 += omm +omx
        elif(self.mod == 3):
            lamdmde = LambdaGeneral(params_info = self.dmde_info)
            omm = self.omegac/a**3
            Ez2 += omm
            Ez2 = Ez2/(1-lamdmde.Omega_x(a))
        else:
            print('Error: no such model. [Hofz]')
            hz = self.H0
            sys.exit(0)
        hz = np.sqrt(Ez2)*self.H0
        return hz

    cpdef double dtauda(self, a):
        '''
        dtauda(a)
        ==========================
        dtauda = c/(a^2 H(z))
        '''
        cdef double dtau, z
        z = 1.0/a -1
        dtau = const['c']/1000/a**2/self.Hofz(z)
        return dtau

    cpdef double rofchi(self, double chi):
        '''
        rofchi(chi)
        ====================
        '''
        cdef double rr
        if(self.flat < 0):
            rr = np.sin(chi)
        elif(self.flat > 0):
            rr = np.sinh(chi)
        else:
            rr = chi
        return rr

    cpdef double DeltaTime(self, double a1, double a2):
        ''' 
        DeltaTime(a1,a2)
        =================================
        param a1: start a
        param a2: end a
        '''
        cdef double dt, tol
        
        dt, tol = integrate.quad(self.dtauda, a1, a2)
        return dt

    cpdef double ComovingRadialDistance(self, double z):
        '''
        ComovingRadialDistance(z)
        ==================================
        '''
        cdef double dC, a

        a = 1.0/(1+z)
        dC = self.DeltaTime(a, 1.0)
        return dC

    cpdef double AngularDiameterDistance(self, double z):
        '''
        AngularDiameterDistance(z)
        ==================================
        '''
        cdef double dA, dC

        dC = self.ComovingRadialDistance(z)
        dA = self.r/(1+z) * self.rofchi(dC/self.r)
        return dA

    cpdef double LuminosityDistance(self, double z):
        '''
        LuminosityDistance(z)
        ===============================
        '''
        cdef double dL
        dL = self.AngularDiameterDistance(z)*(1+z)**2
        return dL

    cpdef double dmz(self, z):
        '''
        dmz(z)
        ===============================
        The comoving angular diameter distance in [Mpc]
        '''
        cdef double dm, chi
        dm = self.AngularDiameterDistance(z)*(1+z)
        return dm

    cpdef double bao_d_v(self, z):
        '''
        bao_d_v(z)
        ===========================
        The bao effectient distance for volume:
        ba0_d_v = [(1+z)^2*dA**2*c*z/H]^{1/3} = [dM**2*c*z/H]^{1/3}
        '''
        cdef double dv3, dv
        dv3 = self.dmz(z)**2 * const['c']/1000 *z /self.Hofz(z)
        dv = dv3**(1.0/3.0)
        return dv

    cpdef double zstar(self):
        '''
        zstar()
        ==============================
        From Hu & Sugiyama
        '''
        cdef double zs, omdmh2

        omdmh2 = self.omch2+self.omnuh2
        zs = 1048*(1+0.00124*self.ombh2**(-0.738))*(1+ \
                (0.0783*self.ombh2**(-0.238)/(1+39.5*self.ombh2**0.763)) * \
                (omdmh2+self.ombh2)**(0.560/(1+21.1*self.ombh2**1.81)))
        return zs

    cpdef double redshift_star(self):
        '''
        redshift_star()
        =========================
        this is z_star
        '''
        cdef double ombh2, ommh2
        cdef double zstar

        ombh2 = self.omegab*(self.H0/100.0)**2
        ommh2 = (self.omegac+self.omegan)*(self.H0/100.0)**2
        #From Hu & Sugiyama
        zstar =  1048*(1+0.00124*ombh2**(-0.738))*(1+
                (0.0783*ombh2**(-0.238)/(1+39.5*ombh2**0.763)) 
                *(ommh2)**(0.560/(1+21.1*ombh2**1.81)))
        return zstar

    cpdef double redshift_drag(self):
        '''
        redshift_drag()
        ===========================
        '''
        cdef double zdrag, ombh2, ommh2, b1, b2

        ombh2 = self.omegab*(self.H0/100)**2
        ommh2 = (self.omegac +self.omegan)*(self.H0/100)**2
        b1 = 0.313*ommh2**(-0.419) *(1.0+0.607*ommh2**0.674)
        b2 = 0.238*ommh2**0.223
        zdrag = 1291.0*ommh2**0.251*(1.0+b1*ombh2**b2)/(1.0+0.6590*ommh2**0.828)

        return zdrag

    cpdef double dsound_da(self, double a):
        '''
        dsound_da(a)
        ===========================
        approximate form used e.g. by CosmoMC for theta
        rs = \int_0^t cs(t) dt/a = \int_0^a cs(a)/(a^2 H) da
        cs2 = c^2/3 [ 1+ 3rho_b/4rho_g]^{-1}
        drsda = cs2/(a^2 H^2)
        param a: scale factor
        '''
        cpdef double dsda, R, cs

        R=3.0e4*a*self.omegab*(self.H0/100)**2
        # R = 3*grhob*a / (4*grhog) # above is mostly within 0.2% and 
        # used for previous consistency
        cs = 1.0/np.sqrt(3*(1+R))
        dsda = self.dtauda(a)*cs
        return dsda

    cpdef double sound_horizon(self, double a0, double a1):
        '''
        sound_horizon(a0,a1)
        ==============================
        param a0: integrate lower limit
        param a1: integrate upper limit
        '''
        cdef double rs, tol

        rs, atol = integrate.quad(self.dsound_da, a0, a1)
        return rs

    cpdef double theta_star(self):
        '''
        theta_star()
        ============================
        '''
        cdef double zstar, astar, rs, DA, tt
        
        zstar = self.redshift_star()
        astar = 1.0/(1+zstar)
        rs = self.sound_horizon(1e-8,astar)
        DA = self.AngularDiameterDistance(zstar)/astar
        tt = rs/DA
        return tt
    
    cpdef double theta_drag(self):
        '''
        theta_drag
        ============================
        '''
        cdef double zd, ad, rs, DA, tt
        
        zd = self.redshift_drag()
        ad = 1.0/(1+zd)
        rs = self.sound_horizon(1e-8,ad)
        DA = self.AngularDiameterDistance(zd)/ad
        tt = rs/DA
        return tt
