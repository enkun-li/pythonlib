#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: modules.pyx
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-04-09 12:09:47
#==================================

import sys
if(not sys.path[0][-10:-1] == 'pythonlib'):
    sys.path.insert(0, '/home/ekli/myworks/pythonlib/')
from constants import constants as const
from equations import LambdaGeneral
import numpy as np
import yaml
from scipy import integrate

#=======================================================
cdef class ModelParams:
    '''
    ModelParams(paramsfile=None,paramsinfo=None,feedback=0)
    ===============================================================
    You need to support a info file or a dict of parameters infos:
    =>>> paramsfile=yourinfofile or paramsinfo=yourinfodict
    =>>> feedback>0 to print some information
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    The info file or dict should include the following informations:
    param use_physical: True
        param h: H0/100 => Hubble constant
        param ombh2: \Omega_b h^2 => baryon
        param omch2: \Omega_c h^2 => cold dark matter
        param omnuh2: \Omega_\nu h^2 => neutrinos
        param omegak: \Omega_k => spatial curvature
    param use_physical: False
        param H0: H0 => Hubble constant
        param omegab: \Omega_b => baryon
        param omegac: \Omega_c => CDM
        param omegan: \Omega_\nu => neutrinos
        param omegav: \Omega_x => dark energy
    param Tcmb: CMB temperature
    param omegam: present dimensionless energy density of matter (include dark matter, baryon)
    param omegar: present dimensionless energy density of radation
    param w_lam: p/rho for the dark energy (assumed constant -1)
    param wa_ppf: for wCDM model
    param cs2_lam: sound speed of dark energy
    param nps(:): array like, model parameters
    param mod: different dark energy models or 
        0: cosmological constant model (LCDM)
        1: constant EOS model (wCDM)
        2: CPL model (wwaCDM)
        3: early dark energy model
    '''
    cdef double h, ombh2, omch2, omnuh2, Tcmb
    cdef double w_lam, wa_ppf, cs2_lam
    cdef double H0, omegab, omegac, omegav, omegan, omegag
    cdef double omegam, omegar, omegak
    cdef double omeganunomass, omeganumassive
    cdef double[:] nps # new model parameters
    cdef int mod
    cdef double r, curv
    cdef int flat, Ksign
    
    def __cinit__(self, paramsinfo=None, paramsfile=None, int feedback=0):
        cdef int i
        cdef double x, a, z

        # if feedback > 0 
        if(feedback>0):
            print("Initial model parameters")
            if(paramsfile is not None):
                print("parameters file is: %s"%paramsfile)
            elif(paramsinfo is not None):
                print("using parameters info dict")
            else:
                print("There is no info file or info dict. [ModelParams]")
                sys.exit(0)

        # read in parameters info
        if(paramsfile is not None):
            ff = open(paramsfile)
            info = yaml.load(ff)
            ff.close
        elif(paramsinfo is not None):
            info = paramsinfo
        else:
            print("There is no info file or info dict. [ModelParams]")
            sys.exit(0)

        # determine the dark energy model
        self.mod = info['dark_energy_model']
        # read in model parameters
        use_physical = info['use_physical']
        if(use_physical):
            self.h = info['h']
            self.H0 = self.h*100
            self.ombh2 = info['ombh2']
            self.omch2 = info['omch2']
            self.omnuh2 = info['omnuh2']
            self.omegak = info['omegak']
            self.omegab = self.ombh2/(self.h)**2
            self.omegac = self.omch2/(self.h)**2
            self.omegan = self.omnuh2/(self.h)**2
            self.omegav = 1-self.omegak-self.omegab-self.omegac-self.omegan
        else:
            self.H0 = info['H0']
            self.h = self.H0/100
            self.omegab = info['omegab']
            self.omegac = info['omegac']
            self.omegan = info['omegan']
            self.omegav = info['omegav']
            self.omegak = 1-self.omegab-self.omegac-self.omegav-self.omegan
            self.ombh2 = self.omegab*(self.H0/100)**2
            self.omch2 = self.omegac*(self.H0/100)**2
            self.omnuh2 = self.omegan*(self.H0/100)**2
        self.Tcmb=info['Tcmb']
        self.w_lam=info['w_lam']
        self.wa_ppf=info['wa_ppf']
        self.cs2_lam=info['cs2_lam']

        self.nps = np.array(info['nps'], dtype=np.float64)

        self.setup_moduleparams()
        #====================================================
        if(feedback>0):
            print("="*40)
            print("%25s = %12.6f"%('Om_b h^2', self.ombh2))
            print("%25s = %12.6f"%('Om_c h^2', self.omch2))
            print("%25s = %12.6f"%('Om_nu h^2', self.omnuh2))
            print("%25s = %12.6f"%('Om_x', self.omegav))
            print("%25s = %12.6f"%('Om_k', self.omegak))
            print("%25s = %12.6f"%('Om_g', self.omegag))
            print("%25s = %12.6f"%('Om_nu,nomass', self.omeganunomass))
            print("%25s = %12.6f"%('Om_r', self.omegar))
            print("%25s = %12.6f"%('Om_m (1-Om_k-Om_x)', 1-self.omegak-self.omegav))
            print("%25s = %12.6f"%('Om_m (Om_b+Om_c+Om_nu,massive)', self.omegab+self.omegac+self.omegan))
            print("~"*40)
            z= self.zstar()
            print("%25s = %12.6f"%('z_star', z))
            print("%25s = %12.6f"%('r_star', self.sound_horizon(z)))
            print("%25s = %12.6f"%('DM_star [Gpc]', self.AngularDiameterDistance( z )*(1+z) /1000 ))
            print("%25s = %12.6f"%('100 theta', 100*self.cosmoTheta()) )
            demodel = ['LCDM', 'wCDM', 'wwaCDM', 'EDECDM']
            print("%25s : %12s"%('cosmological model', demodel[self.mod-1]))
            flat = ['closed', 'flat', 'open']
            print("%25s : %12s"%('The universe is', flat[self.flat+1]))
            for i in range(2):
                x = -7.6+7.6*i
                a = np.exp(x)
                z = np.exp(-x)-1
                print("%15s%9.3f) = %12.6f"%('H(z=',z, self.Hofz(z)*const['c']/1000) )
                print("%15s%9.3f) = %12.6f"%('d_C(z=',z, self.ComovingRadialDistance(z)) )
                print("%15s%9.3f) = %12.6f"%('d_A(z=',z, self.AngularDiameterDistance(z)) )
                print("%15s%9.3f) = %12.6f"%('dcsda(z=',z, self.dsound_da(a)) )
            print("~"*40)
        return 

    def setup_moduleparams(self):
        # flat, closed, open universe
        if(self.omegak > 1e-8):
            #open
            self.curv = -self.omegak/((const['c']/1000)/self.H0)**2
            self.r = 1.0/np.sqrt(np.abs(self.curv))
            self.Ksign = -1
            self.flat = 1
        elif(self.omegak < -1e-8):
            #closed
            self.curv = -self.omegak/((const['c']/1000)/self.H0)**2
            self.r = 1.0/np.sqrt(np.abs(self.curv))
            self.Ksign = 1
            self.flat = -1
        #if(np.abs(self.omegak) < 1e-8 ):
        else:
            #flat
            self.curv = 0.0
            self.r = 1.0
            self.Ksign = 0
            self.flat = 0
        
        #radation energy density
        self.omegag =  2.38e-5*(self.Tcmb/2.7)**4/self.h**2
        self.omeganunomass = 1.62e-5*(self.Tcmb/2.7)**4/self.h**2
        self.omegar = 4.00e-5*(self.Tcmb/2.7)**4/self.h**2
        self.omegam = self.omegab+self.omegac
        
        return

    cpdef double Hofz(self, double z):
        '''
        Hofz(z)
        ===========================
        Hubble parameter in redshift z [H/c light speed in km/s]
        param z: redshift
        '''
        cdef double a, hz
        a = 1.0/(1+z)
        hz = 1.0/a**2/self.dtauda(a)
        return hz

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
        R=3.0e4*a*self.ombh2
        # R = 3*grhob*a / (4*grhog) # above is mostly within 0.2% and used for previous consistency
        cs = 1.0/np.sqrt(3*(1+R))
        dsda = self.dtauda(a)*cs
        return dsda

    cpdef double sound_horizon(self, double z, double ze=1e8):
        '''
        sound_horizon(z)
        ==============================
        param z: end redshift
        param ze: start scale factor
        '''
        cdef double rs, a, ae, tol

        a = 1/(1.0+z)
        ae = 1/(1.0+ze)

        rs, tol = integrate.quad(self.dsound_da,ae,a)
        return rs

    cpdef double cosmoTheta(self):
        '''
        cosmoTheta()
        ============================
        '''
        cdef double zs, astar, rs, DA, tt
        
        zs = self.zstar()
        astar = 1.0/(1+zs)
        rs = self.sound_horizon(zs)
        DA = self.AngularDiameterDistance(zs)/astar
        tt = rs/DA
        return tt
    
    #========================================================
    cpdef double dtauda(self, double a):
        '''
        dtauda(...)
        ========================
        General cosmological model, for dark matter and dark energy, etc.
                '''
        lam = LambdaGeneral(H0=self.H0,omm=self.omegam,omr=self.omegar,
                omk=self.omegak,omx=self.omegav,nps=self.nps,
                w_lam=self.w_lam,wa_ppf=self.wa_ppf,cs2_lam=self.cs2_lam,
                mod=self.mod)
        dtda = lam.dtauda(a)
        return dtda
