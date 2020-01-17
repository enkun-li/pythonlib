##/usr/bin/env python
#-*- coding: utf-8 -*-  
#==================================
# File Name: neutrinos.py
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-04-06 10:19:23
#==================================

#========================================
# Massive neutrinos
#
import numpy as np

cdef class MassiveNu:
    '''
    MassiveNu()
    =================================================
    const = int q^3 F(q) dq = 7/120*pi^4
    self.const = 7.0/120*np.pi**4
    self.const2 = 5.0/7/np.pi**2
    self.zeta3 = 1.2020569031595942853997
    self.zeta5 = 1.0369277551433699263313
    self.zeta7 = 1.0083492773819228268397
    # zeta3*3/2/pi^2*4/11*((k_B*COBE_CMBTemp/hbar/c)^3* 8*pi*G/3/(100*km/s/
    # megaparsec)^2/(c^2/eV)
    self.neutrino_mass_fac=94.07
    self.nrhopn=2000
    self.am_min = 0.01 #0.02
    #smallest a*m_nu to integrate distribution function rather than using series
    self.am_max = 600.0
    #max a*m_nu to integrate
    self.am_minp=self.am_min*1.1
    self.am_maxp=self.am_max*0.9
    
    #Sample for massive neutrino momentum
    #These settings appear to be OK for P_k accuate at 1e-3 level
    self.nqmax0=80 #maximum array size of q momentum samples
    self.nu_q=np.zeros(self.nqmax0)
    self.nu_int_kernel=np.zeros(self.nqmax0)
    nqmax: actual number of q modes evolves
    #public const,Nu_Init,Nu_background, Nu_rho, Nu_drho,  nqmax0, nqmax, &
    #    nu_int_kernel, nu_q, sum_mnu_for_m1, neutrino_mass_fac
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    double[:] r1,p1,dr1,dp1,ddr1,
    '''
    cdef double const, const2, zeta3, zeta5, zeta7
    cdef double neutrino_mass_fac, am_min, am_max, am_minp, am_maxp
    cdef int nrhopn, nqmax0, nqmax
    cdef double dlnam
    cdef double hierarchy_parameter, summu
    cdef double[:] r1,p1,dr1,dp1,ddr1, nu_q, nu_int_kernel
    cdef int Nu_mass_eigenstates

    def __cinit__(self, int Nu_mass_eigenstates):
        #const = int q^3 F(q) dq = 7/120*pi^4
        self.const = 7.0/120*np.pi**4
        self.const2 = 5.0/7/np.pi**2
        self.zeta3 = 1.2020569031595942853997
        self.zeta5 = 1.0369277551433699263313
        self.zeta7 = 1.0083492773819228268397
        # zeta3*3/2/pi^2*4/11*((k_B*COBE_CMBTemp/hbar/c)^3* 8*pi*G/3/(100*km/s/
        # megaparsec)^2/(c^2/eV)
        self.neutrino_mass_fac=94.07
        self.nrhopn=2000
        self.am_min = 0.01 #0.02
        #smallest a*m_nu to integrate distribution function rather than using series
        self.am_max = 600.0
        #max a*m_nu to integrate
        self.am_minp=self.am_min*1.1
        self.am_maxp=self.am_max*0.9
        
        #Sample for massive neutrino momentum
        #These settings appear to be OK for P_k accuate at 1e-3 level
        self.nqmax0=80 #maximum array size of q momentum samples
        self.nu_q=np.zeros(self.nqmax0)
        self.nu_int_kernel=np.zeros(self.nqmax0)

        #Number of Nu_mass_eigenstates
        self.Nu_mass_eigenstates = Nu_mass_eigenstates
        return

    def sum_mnu_for_m1(self, double Delta, double D31sq, double D21sq):
        '''
        sum_mnu_for_m1(Delta, D31sq, D21sq)
        ======================================
        Delta: hierarchy_parameter
        D31sq: \Delta m_{31}^2 = m_3^2-m_1^2
        D21sq: \Delta m_{21}^2 = m_2^2-m_1^2
        '''
        cdef double m1,m2,m3
        m1 = (1-Delta)/2*np.sqrt(D31sq/Delta)
        m2 = np.sqrt((1-Delta)**2/4*D31sq/Delta+D21sq)
        m3 = (1+Delta)/2*np.sqrt(D31sq/Delta)
        self.summnu = m1+m2+m3
        self.hierarchy_parameter = Delta
        return m1,m2,m3

   # def Nu_Init(self):
   #     '''
   #     Nu_Init()
   #     ===========================================
   #     #  Initialize interpolation tables for massive neutrinos.
   #     #  Use cubic splines interpolation of log rhonu and pnu vs. log a*m.
   #     '''
   #     cdef int i
   #     cdef double dq,dlfdlq, q, am, rhonu,pnu
   #     cdef double spline_data(nrhopn)
   #     #  nu_masses=m_nu(i)*c**2/(k_B*T_nu0).
   #     #  Get number density n of neutrinos from
   #     #  rho_massless/n = int q^3/(1+e^q) / int q^2/(1+e^q)=7/180 pi^4/Zeta(3)
   #     #  then m = Omega_nu/N_nu rho_crit /n
   #     #  Error due to velocity < 1e-5
   #     
   #     for i from 0<=i<self.Nu_mass_eigenstates:
   #         nu_masses(i)=const/(1.5d0*zeta3)*grhom/grhor*CP%omegan*CP%Nu_mass_fractions(i) &
   #         /CP%Nu_mass_degeneracies(i)
   # end do

   # if (allocated(r1)) return
   # allocate(r1(nrhopn),p1(nrhopn),dr1(nrhopn),dp1(nrhopn),ddr1(nrhopn))


   # nqmax=3
   # if (AccuracyBoost >1) nqmax=4
   # if (AccuracyBoost >2) nqmax=5
   # if (AccuracyBoost >3) nqmax=nint(AccuracyBoost*10)
   # #note this may well be worse than the 5 optimized points

   # if (nqmax > nqmax0) call MpiStop('Nu_Init: qmax > nqmax0')

   # #We evolve evolve 4F_l/dlfdlq(i), so kernel includes dlfdlnq factor
   # #Integration scheme gets (Fermi-Dirac thing)*q^n exact,for n=-4, -2..2
   # #see CAMB notes
   # if (nqmax==3) then
   #     #Accurate at 2e-4 level
   #     nu_q(1:3) = (/0.913201, 3.37517, 7.79184/)
   #     nu_int_kernel(1:3) = (/0.0687359, 3.31435, 2.29911/)
   # else if (nqmax==4) then
   #     #This seems to be very accurate (limited by other numerics)
   #     nu_q(1:4) = (/0.7, 2.62814, 5.90428, 12.0/)
   #     nu_int_kernel(1:4) = (/0.0200251, 1.84539, 3.52736, 0.289427/)
   # else if (nqmax==5) then
   #     #exact for n=-4,-2..3
   #     #This seems to be very accurate (limited by other numerics)
   #     nu_q(1:5) = (/0.583165, 2.0, 4.0, 7.26582, 13.0/)
   #     nu_int_kernel(1:5) = (/0.0081201, 0.689407, 2.8063, 2.05156, 0.126817/)
   # else
   #     dq = (12 + nqmax/5)/real(nqmax)
   #     do i=1,nqmax
   #         q=(i-0.5d0)*dq
   #         nu_q(i) = q
   #         dlfdlq=-q/(1.+exp(-q))
   #         nu_int_kernel(i)=dq*q**3/(exp(q)+1.) * (-0.25*dlfdlq) #now evolve 4F_l/dlfdlq(i)
   #     end do
   # end if
   # nu_int_kernel=nu_int_kernel/const

   # dlnam=-(log(am_min/am_max))/(nrhopn-1)


   # #$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC) &
   # #$OMP & PRIVATE(am, rhonu,pnu)
   # do i=1,nrhopn
   #     am=am_min*exp((i-1)*dlnam)
   #     call nuRhoPres(am,rhonu,pnu)
   #     r1(i)=log(rhonu)
   #     p1(i)=log(pnu)
   # end do
   # #$OMP END PARALLEL DO


   # call splini(spline_data,nrhopn)
   # call splder(r1,dr1,nrhopn,spline_data)
   # call splder(p1,dp1,nrhopn,spline_data)
   # call splder(dr1,ddr1,nrhopn,spline_data)


   # end subroutine Nu_init

   # #cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   # subroutine nuRhoPres(am,rhonu,pnu)
   # #  Compute the density and pressure of one eigenstate of massive neutrinos,
   # #  in units of the mean density of one flavor of massless neutrinos.

   # real(dl),  parameter :: qmax=30.
   # integer, parameter :: nq=100
   # real(dl) dum1(nq+1),dum2(nq+1)
   # real(dl), intent(in) :: am
   # real(dl), intent(out) ::  rhonu,pnu
   # integer i
   # real(dl) q,aq,v,aqdn,adq


   # #  q is the comoving momentum in units of k_B*T_nu0/c.
   # #  Integrate up to qmax and then use asymptotic expansion for remainder.
   # adq=qmax/nq
   # dum1(1)=0.
   # dum2(1)=0.
   # do  i=1,nq
   #     q=i*adq
   #     aq=am/q
   #     v=1./sqrt(1.+aq*aq)
   #     aqdn=adq*q*q*q/(exp(q)+1.)
   #     dum1(i+1)=aqdn/v
   #     dum2(i+1)=aqdn*v
   # end do
   # call splint(dum1,rhonu,nq+1)
   # call splint(dum2,pnu,nq+1)
   # #  Apply asymptotic corrrection for q>qmax and normalize by relativistic
   # #  energy density.
   # rhonu=(rhonu+dum1(nq+1)/adq)/const
   # pnu=(pnu+dum2(nq+1)/adq)/const/3.

   # end subroutine nuRhoPres

   # #cccccccccccccccccccccccccccccccccccccccccc
   # subroutine Nu_background(am,rhonu,pnu)
   # use precision
   # use ModelParams
   # real(dl), intent(in) :: am
   # real(dl), intent(out) :: rhonu, pnu

   # #  Compute massive neutrino density and pressure in units of the mean
   # #  density of one eigenstate of massless neutrinos.  Use cubic splines to
   # #  interpolate from a table.

   # real(dl) d
   # integer i

   # if (am <= am_minp) then
   #     rhonu=1. + const2*am**2
   #     pnu=(2-rhonu)/3.
   #     return
   # else if (am >= am_maxp) then
   #     rhonu = 3/(2*const)*(zeta3*am + (15*zeta5)/2/am)
   #     pnu = 900./120./const*(zeta5-63./4*Zeta7/am**2)/am
   #     return
   # end if


   # d=log(am/am_min)/dlnam+1.
   # i=int(d)
   # d=d-i

   # #  Cubic spline interpolation.
   # rhonu=r1(i)+d*(dr1(i)+d*(3.*(r1(i+1)-r1(i))-2.*dr1(i) &
   #     -dr1(i+1)+d*(dr1(i)+dr1(i+1)+2.*(r1(i)-r1(i+1)))))
   # pnu=p1(i)+d*(dp1(i)+d*(3.*(p1(i+1)-p1(i))-2.*dp1(i) &
   #     -dp1(i+1)+d*(dp1(i)+dp1(i+1)+2.*(p1(i)-p1(i+1)))))
   # rhonu=exp(rhonu)
   # pnu=exp(pnu)

   # end subroutine Nu_background

   # #cccccccccccccccccccccccccccccccccccccccccc
   # subroutine Nu_rho(am,rhonu)
   # use precision
   # use ModelParams
   # real(dl), intent(in) :: am
   # real(dl), intent(out) :: rhonu

   # #  Compute massive neutrino density in units of the mean
   # #  density of one eigenstate of massless neutrinos.  Use cubic splines to
   # #  interpolate from a table.

   # real(dl) d
   # integer i

   # if (am <= am_minp) then
   #     rhonu=1. + const2*am**2
   #     return
   # else if (am >= am_maxp) then
   #     rhonu = 3/(2*const)*(zeta3*am + (15*zeta5)/2/am)
   #     return
   # end if

   # d=log(am/am_min)/dlnam+1.
   # i=int(d)
   # d=d-i

   # #  Cubic spline interpolation.
   # rhonu=r1(i)+d*(dr1(i)+d*(3.*(r1(i+1)-r1(i))-2.*dr1(i) &
   #     -dr1(i+1)+d*(dr1(i)+dr1(i+1)+2.*(r1(i)-r1(i+1)))))
   # rhonu=exp(rhonu)
   # end subroutine Nu_rho

   # #ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   # function Nu_drho(am,adotoa,rhonu) result (rhonudot)
   # use precision
   # use ModelParams

   # #  Compute the time derivative of the mean density in massive neutrinos
   # #  and the shear perturbation.
   # real(dl) adotoa,rhonu,rhonudot
   # real(dl) d
   # real(dl), intent(IN) :: am
   # integer i

   # if (am< am_minp) then
   #     rhonudot = 2*const2*am**2*adotoa
   # else if (am>am_maxp) then
   #     rhonudot = 3/(2*const)*(zeta3*am - (15*zeta5)/2/am)*adotoa
   # else
   #     d=log(am/am_min)/dlnam+1.
   #     i=int(d)
   #     d=d-i
   #     #  Cubic spline interpolation for rhonudot.
   #     rhonudot=dr1(i)+d*(ddr1(i)+d*(3.*(dr1(i+1)-dr1(i)) &
   #         -2.*ddr1(i)-ddr1(i+1)+d*(ddr1(i)+ddr1(i+1) &
   #         +2.*(dr1(i)-dr1(i+1)))))

   #     rhonudot=rhonu*adotoa*rhonudot/dlnam
   # end if

   # end function Nu_drho

   # end module MassiveNu

   # # wrapper function to avoid cirular module references
   # subroutine init_massive_nu(has_massive_nu)
   # use MassiveNu
   # use ModelParams
   # implicit none
   # logical, intent(IN) :: has_massive_nu

   # if (has_massive_nu) then
   #     call Nu_Init
   # else
   #     nu_masses = 0
   # end if
   # end subroutine in it_massive_nu


#===============================
if __name__ == "__main__":
    mnu = MassiveNu()
    print(mnu.const) 

