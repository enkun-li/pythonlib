#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: cosmography.pyx
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-02-24 08:45:39
#==================================

cimport cython
import sys
# This is the python lib of yourself, so change it with your own dir.
if(not sys.path[0][-10:-1] == 'pythonlib'):
    sys.path.insert(0, '/home/ekli/myworks/pythonlib/')
from darksector import darksector
import numpy as np
import scipy.optimize as opt
from scipy import integrate
from scipy import interpolate

cdef extern from "math.h":
    double sqrt(double xx)
    double exp(double xx)
    double sinh(double xx)
    double cosh(double xx)
    double sin(double xx)
    double cos(double xx)
    double fabs(double xx)
    double log10(double xx)
    double log(double xx)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cpdef double cubic_spline(double x, double[:] xx, double[:] yy,double[:] yy2,
        int order=0):
    '''
    cubic_spline(x, xx, yy, yy2)
    =====================================
    From array {x,f(x),f''(x)} to generate f(z)
    param x: the aim position, steep must be equal!
    param xx: x axis of data
    param yy: f(x) of data
    param yy2: second derivation of f(x) => f''(x)
    param order: 0 for f(x)
                 1 for f'(x)
                 2 for f''(x)
    '''
    cdef int n, j
    cdef double h, A, B
    n = xx.shape[0]
    h = (xx[-1]-xx[0])/(n-1)
    j = int((x-xx[0])/h)
    B = (x-xx[j])/h
    A = (xx[j+1]-x)/h
    if(order == 0):
        fx = A*yy[j] +B*yy[j+1] +h**2/6*( (A**3-A)*yy2[j] +(B**3-B)*yy2[j+1] )
    elif(order==1):
        fx = (yy[j+1]-yy[j])/h -(3*A**2-1)/6*h*yy2[j] +(3*B**2-1)/6*h*yy2[j+1]
    elif(order==2):
        fx = A*yy2[j] +B*yy2[j+1]
    else:
        print('Can not calculate %d order.'%order)
    #print(n, h, j)
    return fx

cpdef double cubic_spline_1d(double x,double[:] xx, double[:] yy, double[:] dyy):
    '''
    cubic_spline_1d(x, xx, yy, dyy)
    =====================================
    From array {x,f(x),f''(x)} to generate f(z)
    param x: the aim position, steep must be equal!
    param xx: x axis of data
    param yy: f(x) of data
    param dyy: first derivation of f(x) => f'(x)
    '''
    cdef int n, j
    cdef double h, dx, A, B
    n = xx.shape[0]
    h = (xx[-1]-xx[0])/(n-1)
    j = int((x-xx[0])/h)
    dx = x - xx[j]
    A = (3*yy[j+1] - 3*yy[j] -(dyy[j+1]+2*dyy[j])*h )/h**2
    B = ( (dyy[j+1]+dyy[j])*h + 2*(yy[j] -yy[j+1]) )
    fx = yy[j] +dyy[j]*dx +A *dx**2 +B *dx**3
    return fx

cpdef double[:] spline(double[:] xx, double[:] yy, double yp0=0, double ypn=0):
    '''
    spline(xx,yy,yp0,ypn)
    =============================
    generate the second dervation of array.
    '''
    cdef int n, k, i
    cdef double[:] y2, u
    cdef double sig, p, qn, un

    n = xx.shape[0]
    y2 = np.zeros(n)
    u = np.zeros(n+1)
    if(yp0>1e30):
        y2[0]=0
        u[0] = 0
    else:
        y2[0] = -0.5
        u[0] = (3/(xx[1] -xx[0]))* ((yy[1]-yy[0])/(xx[1]-xx[0]) -yp0)
    for i in range(1,n-1):
        sig = (xx[i]-xx[i-1])/(xx[i+1]-xx[i-1])
        p = sig*y2[i-1] +2
        y2[i] = (sig-1)/p
        u[i] = (6.*((yy[i+1]-yy[i])/(xx[i+1]-xx[i])-(yy[i]-yy[i-1])/
            (xx[i]-xx[i-1]))/(xx[i+1]-xx[i-1])-sig*u[i-1])/p
    if(ypn > 1e30):
        qn =0
        un = 0
    else:
        qn = 0.5
        un = (3/(xx[n-1] -xx[n-2]))*(ypn-(yy[n-1]-yy[n-2])/(xx[n-1]-xx[n-2]))
        y2[n-1] = (un-qn*u[n-2])/(qn*y2[n-2]+1)
    for k in np.arange(n-2,0,-1):
        y2[k] = y2[k]*y2[k+1] +u[k]        
    return y2

#******************************
# LCDM model
#
cdef class cosmolgoy_bk:
    '''
    cosmolgoy_bk(H0,omm,omb,omk,w0,wa,Tcmb,N_eff,nps,MOD,feedback)
    =================================================
    The general LCDM model: (DM & lambda)
    - H(z)^2 = \kappa (\rho_m + \rho_x)
    - rhom = H0^2/kappa * omm*(1+z)^3
    - rhox = H0^2/kappa * (1-omm)
    ======================================
    The extened LCDM model: wCDM
    - H^2/H0^2 = omm*(1+z)^3 + (1-omm)*(1+z)^(3*(1+w))
    ======================================
    The extened LCDM model: w0waCDM
    - a = 1/(1+z)
    - E^2 = omm*a^{-3} + (1-omm)*a^{-3(1+w0+wa)}*e^{-3wa(1-a)}
    ======================================
    Extened LCDM model with spatial curvature:
    - omk = - Kc^2/H0^2/a_0^2
    - Omk = omk*(1+z)^2
    - E^2 = Omm + Omx + Omk
    ======================================
    kappa = 8\piG/3c^2
    =>>> initial the class:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * Input: 
    H0    : hubble constant
    omm   : matter energy density, from this two, we can derive 
            omc = omm - omb, which is the dark matter energy density
    omb   : baryon energy density
    w0=-1 : present value of w_x
    wa=0  : -dw/da at present
    omk=0 : spatial curvature
    Tcmb=2.7255 : temperature of CMB, this imposed that
            \Omega_\gamma h^2 = 2.38e-5(Tcmb/2.7)^4
    N_eff : effective numbers of neutrinos (only massless neutrinos considered)
    nps   : array of new parameters
    mod   : model of dark energy
    feedback: 0 --> no print, 1 --> print some thing
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    mod:=>> [0,1,2,3] -> [lcdm, wcdm, wwacdm, edecdm]
        =>> [4] -> [IDE01]
    '''
    cdef double H0, omb, omm, omx, omk, w0, wa
    cdef double omg, omnu, omr, omc
    cdef double fK
    cdef double c, Tcmb, N_eff
    cdef double[:] nps
    cdef int mod
    cdef double[:] xx, fc, fx, dfc, dfx, d2fc, d2fx # for differential equations to generate DM & DE rho

    def __cinit__(self, double H0, double omm, double omb=0.04,
            double w0=-1.0, double wa=0.0, double omk = 0.0, 
            double Tcmb=2.7255, double N_eff=3.046, nps=None, 
            MOD='lcdm', feedback=0):
        '''
        Input parameter: H0, omm, w, wa, omk, Tcmb
        '''        
        cdef double zstar, rs, zdrag

        self.Tcmb = Tcmb
        self.H0 = H0
        self.omb = omb
        self.omm = omm
        self.omc = omm-omb
        self.omg =  2.38e-5*(Tcmb/2.7)**4/(H0/100.0)**2
        self.omnu = 1.62e-5*(Tcmb/2.7)**4/(H0/100.0)**2 *(N_eff/3.0)
        self.omr = self.omg+self.omnu
        self.omk = omk
        self.omx = 1.0-omm-omk-self.omr
        self.w0 = w0
        self.wa = wa
        self.omk = omk
        self.c = 2.99792458e5 #[km/s]
        self.fK = self.H0*sqrt(abs(self.omk))/self.c if(abs(self.omk) > 1e-8) else 1.0
        if(nps is not None):
            self.nps = nps
        MODS = {'lcdm':0, 'wcdm':1, 'wwacdm': 2, 'edecdm':3, \
                'IDE01': 4, 'IDE02': 5}
        self.mod = MODS[MOD]
        if(self.mod == 5):
            self.__initial_rho_DMDE()            
        else:
            result_rho = np.zeros((3,3))

        if(feedback == 1):
            print('='*40)
            print('Adopting the %s model.'%MOD)
            print('%25s : %12.5f'%('Hubble constant', self.H0))
            print('%25s : %12.5f'%('Baryon energy density', self.omb))
            print('%25s : %12.5f'%('Dark matter', self.omc))
            print('%25s : %12.5f'%('Dark energy', self.omx))
            print('%25s : %12.5f'%('Photon', self.omg))
            print('%25s : %12.5f'%('Radiation', self.omr))
            print('%25s : %12.5f'%('Spatial curavature', self.omk))
            print('~'*40)
            zstar = self.redshift_star()
            print('%25s : %12.5f'%('z_star', zstar))
            print('%25s : %12.5f'%('D_M(z_star) [Gpc]', self.dmz(zstar)/1000))
            #print('%25s : %12.5f'%('D_A(z_star)', self.angulardiameterdistance(zstar)))
            rs = self.sound_horizon(1e-8,1/(1+zstar))
            print('%25s : %12.5f'%('r_s(z_star)', rs))
            rs = self.sound_horizon_exact(1e-8,1/(1+zstar))
            print('%25s : %12.5f'%('r_s(z_star)', rs))
            print('%25s : %12.5f'%('100 theta_star', 100*self.cmbtotheta()))
            print('~'*40)
            zdrag = self.redshift_drag()
            print('%25s : %12.5f'%('z_drag', zdrag))
            print('%25s : %12.5f'%('D_M(z_drag) [Gpc]', self.dmz(zdrag)/1000))
            rs = self.sound_horizon(1e-8,1/(1+zdrag))
            print('%25s : %12.5f'%('r_s(z_drag)', rs))
            rs = self.sound_horizon_exact(1e-8,1/(1+zdrag))
            print('%25s : %12.5f'%('r_s(z_drag)', rs))
        return

    cpdef __initial_rho_DMDE(self):
        '''
        __initial_rho_DM_DE()
        ===================================
        '''
        cdef double[:,:] result_rho
        cdef double yp0, ypn

        dks = darksector(omegac=self.omc, omegax=self.omx,
                nps=self.nps, mod=self.mod)
        result_rho = dks.initial_rho_DM_DE()
        self.xx = np.array(result_rho[0])
        self.fc = np.array(result_rho[1])
        self.fx = np.array(result_rho[2])
        self.dfc = np.array(result_rho[3])
        self.dfx = np.array(result_rho[4])
        #cs = interpolate.CubicSpline(self.xx, self.fc)
        #self.d2fc = np.array(cs(self.xx, 2), dtype=np.float64)
        #cs = interpolate.CubicSpline(self.xx, self.fx)
        #self.d2fx = np.array(cs(self.xx, 2), dtype=np.float64)
        
        yp0 = (self.fc[1]-self.fc[0])/(self.xx[1]-self.xx[0])
        ypn = (self.fc[-1]-self.fc[-2])/(self.xx[-1]-self.xx[-2])
        self.d2fc = spline(self.xx, self.fc, yp0=yp0,ypn=ypn)
        
        yp0 = (self.fx[1]-self.fx[0])/(self.xx[1]-self.xx[0])
        ypn = (self.fx[-1]-self.fx[-2])/(self.xx[-1]-self.xx[-2])
        self.d2fx = spline(self.xx, self.fx, yp0=yp0,ypn=ypn)
        
        return

    cpdef double Omega_x(self, double z):
        '''
        Omega_x(z)
        =========================
        '''
        cdef double omx, ome, omx1, omx2, omx3
        cdef double a
        a = 1.0/(1.0+z)
        omx =0
        if(self.mod==3):
            ome = self.nps[0]
            omx1 = ome*(1-a**(-3*self.w0))
            omx2 = self.omx+self.omm*a**(3*self.w0) +self.omr*a**(3*self.w0-1)
            omx3 = self.omx-omx1
            omx = omx3/omx2+omx1
        elif(self.mod==2 or self.mod==1 or self.mod==0):
            omx1 = self.omm*(1+z)**3 +self.omr*(1+z)**4 +self.omk*(1+z)**2
            omx2 = self.omx*a**(-3*(1+self.w0+self.wa)) *np.exp(-3*self.wa*(1-a))
            omx3 = omx1+omx2
            omx = omx2/omx3
        else:
            print('No such dark energy model, please define it!')
            sys.exit(0)
        return omx

    cpdef double domxda(self, double z):
        '''
        domxda(z)
        ==============================
        '''
        cdef double domx, ome, omx1, omx2, omx3
        cdef double domx1, domx2, domx3
        cdef double a, w0, wa, omm, omr, omk
        w0 = self.w0
        wa = self.wa
        omx = self.omx
        omm = self.omm
        omr = self.omr
        omk = self.omk
        a = 1.0/(1.0+z)
        domx = 0
        if(self.mod==3):
            ome = self.nps[0]
            omx1 = ome*(1-a**(-3*w0))
            omx2 = omx+omm*a**(3*w0) +omr*a**(3*w0-1)
            omx3 = omx-omx1
            domx1 = 3*w0*ome*a**(-3*w0-1)
            domx2 = 3*w0*omm*a**(3*w0-1) +(3*w0-1)*omr*a**(3*w0-2)
            domx3 = -domx1
            domx = domx3/omx2 -omx3*domx2/omx2**2 +domx1
        elif(self.mod==2 or self.mod==1 or self.mod==0):
            omx1 = omm/a**3 +omr/a**4 +omk/a**2
            omx2 = omx*a**(-3*(1+w0+wa)) *np.exp(-3*wa*(1-a))
            omx3 = omx1+omx2
            domx1 = -3*omm/a**4-4*omr/a**5 -2*omk/a**3
            domx2 = -3*(1+w0+wa)*omx*a**(-3*(1+w0+wa)-1) *np.exp(-3*w0*(1-a)) \
                    +omx*a**(-3*(1+w0+wa))*np.exp(-3*wa*(1-a))*3*wa
            domx3 = domx1 +domx2
            domx = domx2/omx3 -omx2*domx3/omx3**2

        return domx

    cpdef double EoS(self, double z):
        '''
        EoS(z)
        =====================
        Equation of state parameter of DE
        '''
        cdef double w
        cdef double omx, ome, omm, omr, omk
        cdef double a
        a = 1.0/(1.0+z)
        omx = self.omx
        omm = self.omm
        omr = self.omr
        w = 0
        if(self.mod==3):
            w = 1.0/3.0*(omr/(omm*a+omr) \
                    -a*self.domxda(z)/self.Omega_x(z)/(1-self.Omega_x(z)))
        elif(self.mod==2 or self.mod==1 or self.mod==0):
            w = self.w0 + self.wa*z/(1.0+z)
        return w

    cpdef double Hofz(self, double z):
        '''
        Hofz(z)
        =======================
        The Hubble parameter in [km/s/Mpc]
        '''
        cdef double a, x, omm, omx, omk, omr, omb, Ez2, hz
        a = 1.0/(1.0+z)
        Ez2 =0
        if(self.mod==3):
            Ez2 = (self.omm/a**3 +self.omr/a**4)/(1.0-self.Omega_x(z))
        elif(self.mod==2 or self.mod==1 or self.mod==0):
            omm = self.omm*(1+z)**3
            omk = self.omk*(1+z)**2
            omr = self.omr*(1+z)**4
            omx = self.omx* a**(-3*(1+self.w0+self.wa)) *exp(-3*self.wa *(1-a))
            Ez2 = omk + omm + omr + omx
        elif(self.mod == 4):
            dks = darksector(omegac=self.omc, omegax=self.omx, \
                    nps=self.nps, mod=self.mod)
            omm = dks.Erhodm(z)
            omx = dks.Erhodx(z)
            omk = self.omk*(1+z)**2
            omb = self.omb*(1+z)**3
            omr = self.omr*(1+z)**4
            Ez2 =  omm + omx + omk + omb + omr
        elif(self.mod == 5):
            x = np.log(a)
            #csc = interpolate.CubicSpline(self.xx, self.fc)
            #csx = interpolate.CubicSpline(self.xx, self.fx)
            
            #csc = interpolate.interp1d(self.xx, self.fc, kind='cubic')
            #csx = interpolate.interp1d(self.xx, self.fx, kind='cubic')
            
            #omm = self.omc * csc(x)
            #omx = self.omx * csx(x)

            #csc = cubic_spline(x, self.xx, self.fc, self.d2fc)
            #csx = cubic_spline(x, self.xx, self.fc, self.d2fx)

            csc = cubic_spline_1d(x, self.xx, self.fc, self.dfc)
            csx = cubic_spline_1d(x, self.xx, self.fx, self.dfx)
           
            omm = self.omc * exp(csc)
            omx = self.omx * exp(csx)
           
            omk = self.omk*(1+z)**2
            omb = self.omb*(1+z)**3
            omr = self.omr*(1+z)**4
            Ez2 =  omm + omx + omk + omb + omr            

        hz = sqrt(Ez2)*self.H0
        return hz 

    cpdef double omega_m(self, double z):
        '''
        omega_m(z)
        =======================
        Dimensionless energy density of DM
        '''
        cdef double omi, omm
        omi = self.omm*(1+z)**3
        omm = omi/self.Hofz(z)**2
        return omm

    #cpdef double omega_x(self, double z):
    #    '''
    #    omega_x(z)
    #    =======================
    #    Dimensionless energy density of DE
    #    '''
    #    cdef double a, omi, omx
    #    a = 1./(1.+z)
    #    omi = self.omx* a**(-3*(1+self.w0+self.wa)) *exp(-3*self.wa *(1-a))
    #    omx = omi/self.Hofz(z)**2
    #    return omx

    cpdef double omega_k(self, double z):
        '''
        omega_k(z)
        ========================
        Dimensionless energy density of spatial curvature
        '''
        cdef double omi, omk
        omi = self.omk*(1+z)**2
        omk = omi/self.Hofz(z)**2
        return omk

    cpdef double dtauda(self, double z):
        '''
        dtauda(z)
        ========================
        dtauda = 1/H(z)
        '''
        cdef double dta
        dta = 1.0/self.Hofz(z)
        return dta

    cpdef double comovingdistance(self, double z, double z0=0.0):
        '''
        comovingdistance(z, z0)
        =================================
        Comoving distance: [Mpc]
        - z : the upper integrate limit
        - z0: the lower integrate limit
        '''
        cdef double dci, dc, err, x, x0
        dci, err = integrate.quad(self.dtauda, z0, z)
        dc = dci*self.c
        return dc

    cpdef double __sinn(self, double chi):
        '''
        sinn(chi)
        =================================
        - sinh(chi) :omk > 0 :K < 0
        - chi       :omk = 0 :K = 0
        - sin(chi)  :omk < 0 :K > 0
        '''
        cdef double sink
        if(self.omk > 1e-8):
            sink = sinh(chi)
        elif(self.omk < -1e-8):
            sink = sin(chi)
        else:
            sink = chi
        return sink

    cpdef double dmz(self, double z):
        '''
        dmz(z)
        ===============================
        The comoving angular diameter distance in [Mpc]
        '''
        cdef double dm, chi
        chi = self.fK*self.comovingdistance(z)
        dm = 1.0/self.fK* self.__sinn(chi)
        return dm

    cpdef double luminositydistance(self, double z):
        '''
        luminositydistance(z)
        ==================================
        The luminosity distance
        '''
        cdef double dl
        dl = (1+z)*self.dmz(z)
        return dl

    cpdef double angulardiameterdistance(self, double z):
        '''
        angulardiameterdistance(z)
        ==================================
        The angular diameter distance
        '''
        cdef double da
        da = self.dmz(z)/(1+z)
        return da

    cpdef double bao_d_v(self, double z):
        '''
        bao_d_v(z)
        ==========================
        The bao effectient distance for volume:
        ba0_d_v = [(1+z)^2*dA**2*c*z/H]^{1/3} = [dM**2*c*z/H]^{1/3}
        '''
        cdef double dv3, dv
        dv3 = self.dmz(z)**2 * self.c *z /self.Hofz(z)
        dv = dv3**(1.0/3.0)
        return dv
 
    cpdef double distancemodulus(self, double z):
        '''
        distancemodulus(z)
        =======================================
        SN Ia: distance modulus:
        - mu = 5 log10(dL) +25
        '''
        cdef double mu
        mu = 5.0*log10(self.luminositydistance(z)) +25
        return mu  

    cpdef double dsound_da_exact(self, double a):
        '''
        dsound_da_exact(a)
        =========================
        '''
        cdef double z, R, cs, dsda
        z = 1.0/a -1
        R = 3*self.omb*a/(4*self.omg)
        cs = self.c/np.sqrt(3*(1+R))
        dsda = self.dtauda(z)*cs/a**2
        return dsda
    
    cpdef double dsound_da(self, double a):
        '''
        dsound_da(a)
        ==========================
        '''
        cdef double z, R, cs, dsda

        z = 1.0/a-1
        R = 3.0e4*a*self.omb*(self.H0/100.0)**2
        # R = 3*grhob*a / (4*grhog) //above is mostly within 0.2% 
        # and used for previous consistency
        cs = self.c/np.sqrt(3*(1+R))
        dsda = self.dtauda(z)*cs/a**2
        return dsda

    cpdef double redshift_star(self):
        '''
        redshift_star()
        =========================
        this is z_star
        '''
        cdef double ombh2, ommh2
        cdef double zstar

        ombh2 = self.omb*(self.H0/100.0)**2
        ommh2 = self.omm*(self.H0/100.0)**2
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

        ombh2 = self.omb*(self.H0/100)**2
        ommh2 = self.omm*(self.H0/100)**2
        b1 = 0.313*ommh2**(-0.419) *(1.0+0.607*ommh2**0.674)
        b2 = 0.238*ommh2**0.223
        zdrag = 1291.0*ommh2**0.251*(1.0+b1*ombh2**b2)/(1.0+0.6590*ommh2**0.828)

        return zdrag

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

    cpdef double sound_horizon_exact(self, double a0, double a1):
        '''
        sound_horizon(a0,a1)
        ==============================
        param a0: integrate lower limit
        param a1: integrate upper limit
        '''
        cdef double rs, tol

        rs, atol = integrate.quad(self.dsound_da_exact, a0, a1)
        return rs

    cpdef double cmbtotheta(self):
        '''
        cmbtotheta()
        ===========================
        '''
        cdef double zstar, astar, rs, DA
        cdef double theta
        
        zstar = self.redshift_star()
        astar = 1/(1.0+zstar)
        rs = self.sound_horizon(1e-8, astar)
        DA = self.angulardiameterdistance(zstar)*(1+zstar)
        theta = rs/DA
        return theta


#*****************************
# Cosmography model:
#
cdef class cosmography_series:
    '''
    This is for cosmography model.
    - H(z) is calculated by Taylor series of a(t) in y
      where y = 1-a = z/(1+z)
    - dA and dL also calculated by series of a in y
    ================================
    =>>> initial class:
    * Input: parameters, flat=True
    - Initialize the cosmography model.
            :parameters: the model parameters, the max order is 4. 
                The order of parameters should be H0, q0, j0, s0, l0, omk 
                (if the sparial curvature is not flat).
            :flat: whether it is falt spatial curvature or not: 
                if False, the last parameter must be Omega_K; 
                if True, do not care it, just include your parameters in parameters.
    '''
    cdef double[:] pars
    cdef int order
    cdef double c

    def __cinit__(self, double[:] parameters, flat=True, out=False):
        cdef int i, num        
        self.order = 0
        num = np.shape(parameters)[0]
        self.pars = np.array([0.0]*6)
        if flat:
            self.order = num-1
            for i from 0<=i<num:
                self.pars[i] = parameters[i]
        elif (not flat):
            self.order = num-2
            self.pars[-1]=parameters[-1]
            for i from 0<=i<num-1:
                self.pars[i] = parameters[i]
        self.c = 2.99792458e5 #[km/s]
        if out:
            print([self.pars[i] for i from 0<=i<6])
        return
    
    def __par2parameter(self):
        H0 = self.pars[0]
        q0 = self.pars[1]
        j0 = self.pars[2]
        s0 = self.pars[3]
        l0 = self.pars[4]
        omk = self.pars[5]
        return H0, q0, j0, s0, l0, omk


    cpdef double Hofz(self, double z):
        '''
        Hofz(z)
        ===========================
        The Hubble parameter. [km/s/Mpc]
            :z: the redshift.
        '''
        cdef double H0, q0, j0, s0, l0, omk
        cdef double y, hh
        cdef double[:] exps
        cdef int i
        H0, q0, j0, s0, l0, omk = self.__par2parameter()
        y = z/(1.+z)
        exps = np.array([0.0]*5)
        exps[0] = 1
        exps[1] = 1+q0
        exps[2] = 1./2.*(2.+j0+2*q0-q0**2)
        exps[3] = 1./6.*(6+j0*(3-4*q0)+6*q0 -3*q0**2+3*q0**3-s0)
        exps[4] = -1./24.*(-24 + 4*j0**2- l0 -24*q0 + 12*q0**2 - 12*q0**3 \
                + 15*q0**4 + j0* (-12 + 16*q0 - 25*q0**2 ) \
                + 4*s0 - 7*q0*s0)
        
        hh = 0.
        for i from 0 <= i < (self.order+1):
            hh += H0*exps[i]*y**i
        return hh

    cpdef double angulardiameterdistance(self, double z):
        '''
        angulardiameterdistance(z)
        =================================
        The Angular Diameter Distance.
            :z: the redshift.
        '''
        cdef double H0, q0, j0, s0, l0, omk
        cdef double y, DA
        cdef double[:] exps
        cdef int i
        H0, q0, j0, s0, l0, omk = self.__par2parameter()
        y = z/(1.+z)
        exps = np.array([0.0]*6)
        exps[1] = 1.
        exps[2] = -0.5*(1+q0)
        exps[3] = 1./6.*(-1.-j0 + q0 + 3*q0**2 + omk)     
        exps[4] = 1./24.*(-2 + j0+ 2*q0 + 10*j0*q0 - 3*q0**2 \
                - 15*q0**3 +s0 + (2 -6*q0)* omk)
        exps[5] = 1./120.*( -6+10*j0**2 -l0 +6*q0-9*q0**2 +15*q0**3 \
                +105*q0**4 + j0*(3-10*q0-105*q0**2) -s0 -15*q0*s0 \
                -5*(-1+2*j0+4*q0-9*q0**2)*omk )
        DA = 0.
        for i from 0 <=i < (self.order+2):
            DA += self.c/H0*exps[i]*y**i
        return DA

    cpdef double luminositydistance(self, double z):
        '''
        luminositydistance(z)
        ================================
        The Luminosity Distance.
            :z: the redshift.
        '''
        cdef double H0, q0, j0, s0, l0, omk
        cdef double y, DL
        cdef double[:] exps
        cdef int i
        H0, q0, j0, s0, l0, omk = self.__par2parameter()
        y = z/(1.+z)
        exps = np.array([0.0]*6)
        exps[1] = 1.
        exps[2] = 0.5*(3-q0)
        exps[3] = 1./6.*(11.-j0 -5* q0 + 3*q0**2 + omk)     
        exps[4] = 1./24.*(50 -26*q0 + 21*q0**2 - 15*q0**3 \
                +j0*(-7+10*q0) +s0 + (10 -6*q0)* omk)
        exps[5] = 1./120.*( 274 +10*j0**2 -l0 -154*q0 +141*q0**2 -135*q0**3 \
                +105*q0**4 + j0*( -47 +90*q0-105*q0**2) +9*s0 -15*q0*s0 \
                -5*(-17 +2*j0 + 16*q0-9*q0**2)*omk )
        DL = 0.
        for i from 0<= i <(self.order+2):
            DL += self.c/H0*exps[i]*y**i
        return DL

    cpdef double dmz(self, double z):
        '''
        dmz(z)
        ============================
        The comoving angular diameter distance.
        Here we use angulardiameterdistance to calc. it.
        dM = dA*(1+z)
        '''
        cdef double dm
        dm = self.angulardiameterdistance(z)*(1+z)
        return dm

    cpdef double bao_d_v(self, double z):
        '''
        bao_d_v(z)
        ==========================
        The bao effectient distance for volume:
        ba0_d_v = [(1+z)^2*dA**2*c*z/H]^{1/3} = [dM**2*c*z/H]^{1/3}
        '''
        cdef double dv3, dv
        dv3 = self.dmz(z)**2 * self.c *z /self.Hofz(z)
        dv = dv3**(1.0/3.0)
        return dv
    
    cpdef double distancemodulus(self, double z):
        '''
        distancemodulus(z)
        =======================================
        SN Ia: distance modulus:
        - mu = 5 log10(dL) +25
        '''
        cdef double mu
        mu = 5.0*log10(self.luminositydistance(z)) +25
        return mu

#**************************************************
cdef class ModelCosmoIny:
    '''
    The cosmolography model series in y.
    =====================================
    y = z/(1+z)
    y => (-oo, 1]
    H0, q0, j0, s0
    '''
    cdef double H0, q0, j0, s0, l0, omk
    cdef double c_light
    cdef int order

    def __cinit__(self, double[:] theta, int order=4):
        cdef double[:] tt
        self.c_light = 2.99792458e5 # [km/s]
        self.order = order
        tt = np.array([67.27,-0.55,1.0,-0.35,3.115, 0.0])
        tt[:order] = theta[:order]
        tt[-1] = theta[-1]
        self.H0, self.q0, self.j0, self.s0, self.l0, self.omk = tt
        return

    cpdef double Hofz(self, double z):
        '''
        Hofz(z)
        ==================
        series in 3 order (z=0)
        '''
        cdef double[:] eps
        cdef int i
        cdef double hz, y

        y = z/(1+z)
        hz = 0.0
        eps = np.zeros(5)

        eps[0] = 1.0
        eps[1] = 1.0+self.q0
        eps[2] = 0.5*(2 + self.j0 +2*self.q0 - self.q0**2)
        eps[3] = 1.0/6.0*(6+6*self.q0 -3*self.q0**2 -3*self.q0**3 \
                +self.j0*(3 -4*self.q0) -self.s0)
        eps[4] = -1./24.*(-24 + 4*self.j0**2- self.l0 -24*self.q0 \
                + 12*self.q0**2 - 12*self.q0**3 + 15*self.q0**4 \
                + self.j0* (-12 + 16*self.q0 - 25*self.q0**2 ) \
                + 4*self.s0 - 7*self.q0*self.s0)

        for i from 0 <= i < self.order:
            hz = hz + eps[i] * y**float(i)
        hz = self.H0*hz
        return hz

    cpdef double comovingdistance(self, double z):
        '''
        comovingdistance(z)
        ============================
        comoving distance in order 4 (z=0)^4
        '''
        cdef double[:] eps
        cdef double dc, y
        cdef int i
        
        y = z/(1+z)
        dc = 0.0
        eps = np.zeros(6)

        eps[0] = 0.0
        eps[1] = 1.0
        eps[2] = 0.5*(1 -self.q0)
        eps[3] = 1.0/6.0*(2 -self.j0 -2*self.q0 +3*self.q0**2)
        eps[4] = 1.0/24.0*(6 -6*self.q0 +9*self.q0**2 -15*self.q0**3 \
                - self.j0*(3 -10*self.q0) +self.s0)
        eps[5] = 1.0/120.0*( 120 + 10*self.j0**2 -self.l0 +60*self.q0**3 
                +105*self.q0**4 -96*(1+self.q0)- 5*self.j0*self.q0*(8+21*self.q0)
                -36*(self.j0-self.q0*(2+3*self.q0))-4*self.s0 -15*self.q0*self.s0
                +8*(-9*self.q0**2-15*self.q0**3+self.j0*(3+10*self.q0) +self.s0))

        for i from 0 <= i <= self.order:
            dc = dc+ eps[i] * y**float(i)
        dc = self.c_light * dc / self.H0
        return dc

    cpdef double luminositydistance(self, double z):
        '''
        luminositydistance(z)
        =========================================
        luminosity distance in order 4
        '''
        cdef double[:] eps
        cdef double dl, y
        cdef int i

        y = z/(1+z)
        dl = 0.0
        eps = np.zeros(6)

        eps[0] = 0.0
        eps[1] = 1.0
        eps[2] = 0.5*(3 -self.q0)
        eps[3] = 1.0/6.0*(11 -self.j0 -5*self.q0 +3*self.q0**2 +self.omk)
        eps[4] = 1.0/24.0*(50 +21*self.q0**2 -15*self.q0**3 \
                +self.j0*(-7 +10*self.q0) +self.s0 \
                +10*self.omk -2*self.q0*(13 +3*self.omk))
        eps[5] = 1./120.*( 274 +10*self.j0**2 -self.l0 -154*self.q0 \
                +141*self.q0**2 -135*self.q0**3 +105*self.q0**4 \
                + self.j0*( -47 +90*self.q0-105*self.q0**2) +9*self.s0 \
                -15*self.q0*self.s0 \
                -5*(-17 +2*self.j0 + 16*self.q0-9*self.q0**2)*self.omk )
        for i from 0 <= i <= self.order:
            dl = dl+ eps[i] * y**float(i)
        dl = self.c_light * dl / self.H0
        return dl

    cpdef double angulardiameterdistance(self, double z):
        '''
        angulardiameterdistance(z)
        ========================================
        angular diameter distance in z
        '''
        cdef double[:] eps
        cdef double da, y
        cdef int i
    
        y = z/(1+z)
        da = 0.0
        eps = np.zeros(6)

        eps[0] = 0.0
        eps[1] = 1.0
        eps[2] = -0.5*(1 +self.q0)
        eps[3] = 1.0/6.0*(-1 -self.j0 +self.q0 +3*self.q0**2 +self.omk)
        eps[4] = 1.0/24.0*(-2 -3*self.q0**2 -15*self.q0**3 \
                + self.j0*(1 +10*self.q0) +self.s0 \
                +2*self.omk +self.q0*(2 -6*self.omk))
        eps[5] = 1./120.*( -6+10*self.j0**2 -self.l0 +6*self.q0-9*self.q0**2 \
                +15*self.q0**3 +105*self.q0**4 + \
                self.j0*(3-10*self.q0-105*self.q0**2) -self.s0 -15*self.q0*self.s0 \
                -5*(-1+2*self.j0+4*self.q0-9*self.q0**2)*self.omk )
        for i from 0 <= i <= self.order:
            da = da+ eps[i] * y**float(i)
        da = self.c_light * da / self.H0
        return da

    cpdef double dmz(self, double z):
        '''
        dmz(z)
        ============================
        The comoving angular diameter distance.
        Here we use angulardiameterdistance to calc. it.
        dM = dA*(1+z)
        '''
        cdef double dm
        dm = self.angulardiameterdistance(z)*(1+z)
        return dm

    cpdef double bao_d_v(self, double z):
        '''
        bao_d_v(z)
        ==========================
        The bao effectient distance for volume:
        ba0_d_v = [(1+z)^2*dA**2*c*z/H]^{1/3} = [dM**2*c*z/H]^{1/3}
        '''
        cdef double dv3, dv
        dv3 = self.dmz(z)**2 * self.c_light *z /self.Hofz(z)
        dv = dv3**(1.0/3.0)
        return dv
    
    cpdef double distancemodulus(self, double z):
        '''
        distancemodulus(z)
        =======================================
        SN Ia: distance modulus:
        - mu = 5 log10(dL) +25
        '''
        cdef double mu
        mu = 5.0*log10(self.luminositydistance(z)) +25
        return mu

#******************************
# Cosmography model:
#
cdef class distancefromdm:
    '''
    All the distance are derivated from the comoving
    angular diameter distance (d_M).
    ===================================
    =>>> initial the class:
    * Input: theta(:)
    -- H0, q0, j0, s0, omk = theta
    '''
    cdef double H0, q0, j0, s0, l0, omk
    cdef double c
    cdef int order
    
    def __cinit__(self,double[:] theta, int order=4):
        cdef double[:] tt 
        tt = np.array([67.27,-0.55,1.0,-0.35,3.115,0.0])
        self.order = order
        tt[-1] = theta[-1]
        tt[:order] = theta[:order]
        self.H0 = tt[0]
        self.q0 = tt[1]
        self.j0 = tt[2]
        self.s0 = tt[3]
        self.l0 = tt[4]
        self.omk = tt[5]
        self.c = 2.99792458e5 #[km/s]
        return
    
    cpdef double dmz(self, double z):
        '''
        comovingdA(z)
        ----------------------
        The comoving angular diameter distance to O(y^4) in [Mpc]
        z: redshift
        '''
        cdef double y, dm
        cdef double[:] eps
        cdef int i
        y = z/(1.+z)
        eps = np.zeros(6)
        eps[1] = 1.0
        eps[2] = 0.5*(1-self.q0)
        eps[3] = 1./6.*(2-self.j0 -2*self.q0 +3.*self.q0**2 +self.omk)
        eps[4] = 1./24.*(6-6*self.q0+9*self.q0**2-15*self.q0**3 \
            +self.j0*(-3.0+10*self.q0) +self.s0 +6*(1-self.q0)*self.omk)
        eps[5] = 1./120.*(24 +10*self.j0**2 -self.l0 -24*self.q0 +36*self.q0**2 
                -60*self.q0**3 +105*self.q0**4 +4*self.s0 -15*self.q0*self.s0 
                -self.j0*(12-40*self.q0+105*self.q0**2)
                +(35 -50*self.q0 +45*self.q0**2-10*self.j0)*self.omk)
        dm = 0.0
        for i from 0<=i<=self.order:
            dm = dm + self.c/self.H0*eps[i]*y**float(i)
        return dm
    
    cpdef double __ddmdy(self, double z):
        cdef double y, ddm
        cdef double[:] eps
        cdef int i
        y = z/(1.+z)
        eps = np.zeros(6)
        eps[1] = 1.0
        eps[2] = 0.5*(1-self.q0)
        eps[3] = 1./6.*(2-self.j0 -2*self.q0 +3.*self.q0**2 +self.omk)
        eps[4] = 1./24.*(6-6*self.q0+9*self.q0**2-15*self.q0**3 \
            +self.j0*(-3.0+10*self.q0) +self.s0 +6*(1-self.q0)*self.omk)
        eps[5] = 1./120.*(24 +10*self.j0**2 -self.l0 -24*self.q0 +36*self.q0**2 
                -60*self.q0**3 +105*self.q0**4 +4*self.s0 -15*self.q0*self.s0 
                -self.j0*(12-40*self.q0+105*self.q0**2)
                +(35 -50*self.q0 +45*self.q0**2-10*self.j0)*self.omk)

        ddm = 0.0
        for i from 0<i<=self.order:
            ddm = ddm + self.c/self.H0*eps[i]*y**float(i-1)*float(i)
        return ddm
        
    cpdef double ddmdz(self, double z):
        '''
        ddm/dz = ddm/dy dy/dz = ddm/dy 1/(1+z)^2
        '''
        cdef double ddm
        ddm = self.__ddmdy(z)/(1.0+z)**2
        return ddm
    
    cpdef double Hofz(self, double z):
        '''
        Hofz(z)
        ----------------------
        Hubble parameters from d_M in [km/s/Mpc]
        '''
        cdef double hz
        hz = self.c/self.ddmdz(z)*sqrt(1+self.H0**2*self.omk/self.c**2*self.dmz(z)**2)
        return hz
    
    cpdef double luminositydistance(self, double z):
        '''
        luminositydistance(z)
        -------------------------------
        d_L = (1+z)d_M
        '''
        cdef double dl
        dl = (1.+z)*self.dmz(z)
        return dl
    
    cpdef double angulardiameterdistance(self, double z):
        '''
        angulardiameterdistance(z)
        -------------------------------
        d_A = d_M/(1+z)
        '''
        cdef double da
        da = self.dmz(z)/(1.+z)
        return da

    cpdef double bao_d_v(self, double z):
        '''
        bao_d_v(z)
        ==========================
        The bao effectient distance for volume:
        ba0_d_v = [(1+z)^2*dA**2*c*z/H]^{1/3} = [dM**2*c*z/H]^{1/3}
        '''
        cdef double dv3, dv
        dv3 = self.dmz(z)**2 * self.c *z /self.Hofz(z)
        dv = dv3**(1.0/3.0)
        return dv

    cpdef double distancemodulus(self, double z):
        '''
        distancemodulus(z)
        =======================================
        SN Ia: distance modulus:
        - mu = 5 log10(dL) +25
        '''
        cdef double mu
        mu = 5.0*log10(self.luminositydistance(z)) +25
        return mu
