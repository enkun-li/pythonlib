#!/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: darksector.pyx
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-04-17 18:05:15
#==================================

cimport cython
import sys
import numpy as np
from scipy import integrate
from scipy import interpolate

#============================
# dark sector
cdef class darksector:
    ''' 
    darksector(nps)
    ==============================
    This is for dark matter and dark energy.
    You can define yourself DM & DE model here, and call it in cosmology_bk.
    param omegac: dark matter energy density
    param omegax: dark energy energy density
    param nps: array like, model parameters
    param mod: for differential models
          =>> 4: IDE01
          =>> 5: IDE02 using differential equation
    '''
    cdef double omegac, omegax
    cdef double[:] nps # for general model parameters
    cdef int mod # for DM/DE model
    cdef double[:] xx, fc, fx # for x, f_c, f_x to CubicSpline, to generate DM/DE energy density

    def __cinit__(self, double omegac=0.26, double omegax=0.7, double[:] nps=np.zeros(5), int mod=4):
        self.omegac = omegac
        self.omegax = omegax
        self.mod = mod
        self.nps = nps
        #if(self.mod==5):
        #    self.initial_rho_DM_DE()
        return

    cpdef double[:] get_rho_DM_DE(self):
        '''
        return result array
        '''
        return self.xx

    cpdef double Erhodm(self, double z):
        '''
        Erhodm(z)
        ===================
        Erhodm = Omega_m a^{-3}
        '''
        cdef double omm, x, xi
        xi = self.nps[1]
        if(self.mod == 4):
            omm = self.omegac * (1+z)**(3-xi)
        else:
            omm = 0
        return omm

    cpdef double Erhodx(self, double z):
        '''
        Erhodx(z)
        ===================
        Erhodx = Omega_x a^{-3(1+w)}
        '''
        cdef double ff, omx, wx, xi
        wx = self.nps[0]
        xi = self.nps[1]
        if(self.mod == 4):
            ff = xi * self.omegac/(xi+3*wx)
            omx = (self.omegax + ff)*(1+z)**(3+3*wx) -ff*(1+z)**(3-xi)
        else:
            omx = 0
        return omx

    cpdef double w_x(self, double z):
        '''
        w_x(x)
        =====================
        equation of state parameter for dark energy.
        '''
        cdef double a, w0, wa, wx
        a = 1/(1+z)
        w0 = self.nps[0]
        wa = self.nps[1]
        wx = w0 +wa*(1-a)
        return wx

    cpdef double f_rhox(self, double z):
        '''
        f_rhox(x)
        ======================
        general formation of DE model, used to initial differential equations.
        '''
        cdef double a, w0, wa, rr
        a = 1/(1+z)
        w0 = self.nps[0]
        wa = self.nps[1]
        rr = a**(-3*(1+w0+wa)) * np.exp(-3*wa*(1-a))
        return rr

    cpdef double interact_xi(self, double fc, double fx):
        '''
        interact_xi(fc,fx)
        ========================
        Q = xi_c rhoc +xi_x rhox
        interactions between DM and DE
        '''
        cdef double xi, xic, xix, xit
        cdef double om, ox, sign, mins

        xi = 0
        if(self.mod==5):
            xic = self.nps[2]
            xix = self.nps[3]
            xi = xic*self.omegac *fc + xix*self.omegax *fx
        elif(self.mod ==6):
            #xit = self.nps[2] *np.min((self.omegax*fx)/(self.omegac*fc), (self.omegax*fx)/(self.omegac*fc))
            om = self.omegac*fc
            ox = self.omegax*fx
            sign = 0 if(om==ox) else (om-ox)/np.abs(om-ox)
            mins = np.min([om/ox, ox/om])
            xit = self.nps[2] *sign *mins
            xi = xit*(self.omegac *fc + self.omegax *fx)
        return xi

    cpdef double[:] derive_array(self, double[:] y, double z):
        '''
        derive_array(y,z)
        =====================
        differential equation of conservation law equations for DM and DE
        param z: = z
        param y: array
        return array
        '''
        cdef double[:] dy
        cdef double fc, fx, flu
        dy = np.zeros(2)
        fc = y[0] #np.exp(y[0])
        fx = y[1] #np.exp(y[1])
        flu = 3*self.interact_xi(fc, fx)
        dy[0] = (3*fc +flu/self.omegac)/(1+z)
        dy[1] = (3*(1+self.w_x(z))*fx -flu/self.omegax)/(1+z)
        return dy

    cpdef double[:,:] initial_rho_DM_DE(self, double xmin=-18.5, double xmax=0.0):
        '''
        initial_rho_DM_DE()
        ==================================
        Using the differential equation to initial array of DM & DE energy densities.
        param xmin: the minimal x
        param xmax: the max x
        '''
        #cdef double[:] xx, y0
        cdef double[:,:] sol, result
        cdef double[:] dy, dyc, dyx
        cdef int i, nrho

        nrho = 1001

        xx = np.linspace(xmin, xmax, nrho)
        zz = np.sort(np.exp(-xx)-1)
        #y0 = np.array([-3*xmin, np.log(self.f_rhox(xmin)) ])
        y0 = np.array([1.0,1.0])
        #sol = integrate.odeint(self.derive_array, y0, xx)
        sol = integrate.odeint(self.derive_array, y0, zz)
        
        dyc = np.zeros(nrho)
        dyx = np.zeros(nrho)
        for i from 0<=i<nrho:
            dy = self.derive_array(sol[i,:], zz[i])
            dyc[i] = -(1+zz[i]) *dy[0]
            dyx[i] = -(1+zz[i]) *dy[1]

        result = np.array([-np.log(1+zz), sol.T[0], sol.T[1], dyc, dyx])
        return result

    #cpdef double[:] derive_array(self, double[:] y, double x):
    #    '''
    #    derive_array(x,y)
    #    =====================
    #    differential equation of conservation law equations for DM and DE
    #    param z: = z
    #    param y: array
    #    return array
    #    '''
    #    cdef double[:] dy
    #    dy = np.zeros(2)
    #    dy[0] = -3*y[0] -3*self.interact_xi(y[0], y[1])/self.omegac
    #    dy[1] = -3*(1+self.w_x(x))*y[1] +3*self.interact_xi(y[0], y[1])/self.omegax
    #    return dy

    #cpdef double[:,:] initial_rho_DM_DE(self, double xmin=-13.4, double xmax=0.4):
    #    '''
    #    initial_rho_DM_DE()
    #    ==================================
    #    Using the differential equation to initial array of DM & DE energy densities.
    #    param xmin: the minimal x
    #    param xmax: the max x
    #    '''
    #    #cdef double[:] xx, y0
    #    cdef double[:,:] sol, result
    #    cdef double[:] dy, dyc, dyx
    #    cdef int i, nrho

    #    nrho = 1001

    #    xx = np.linspace(xmin, xmax, nrho)
    #    y0 = np.array([np.exp(-3*xmin), self.f_rhox(xmin)])
    #    sol = integrate.odeint(self.derive_array, y0, xx)

    #    dyc = np.zeros(nrho)
    #    dyx = np.zeros(nrho)
    #    for i from 0<=i<nrho:
    #        dy = self.derive_array(sol[i,:], xx[i])
    #        dyc[i] = dy[0]
    #        dyx[i] = dy[1]

    #    result = np.array([xx, sol.T[0], sol.T[1], dyc, dyx])
    #    return result
