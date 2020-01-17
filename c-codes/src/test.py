#!/home/ekli/.local/anaconda3/bin/python
#-*- coding: utf-8 -*-  
#==================================
# File Name: test.py
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-06-11 11:26:43
#==================================

import sys
if(not sys.path[0][-10:-1] == 'pythonlib'):
    sys.path.insert(0, '/home/ekli/myworks/pythonlib/')

from my_gapp import *

from gp_likelihood import min_GaussianProcesses as mGP
from general_fsigma import GaussianProcesses as fs8GP

#===================================================
style = {'sty': ['-r','--g', '-.b', ':k','.-r',
                 '-g','--b', '-.k', ':r','.-g',
                 '-b','--k', '-.r', ':g','.-b',
                 '-k','--r', '-.g', ':b','.-k'],
                 }
#===================================================
def readin_dat(datfile):
    rec_dat = np.loadtxt(datfile, unpack=True)
    X = rec_dat[0]
    Y = rec_dat[1]
    err = rec_dat[2]
    cov = np.eye(len(err))*err**2
    return X, Y, cov

def cov_fs8_star(fs8_cij_w, errf):
    cov_w = np.loadtxt(fs8_cij_w, unpack=True)
    n = len(errf)
    cov_star = np.eye(n) * errf**2
    cov_star[9:12,9:12] = cov_w
    return cov_star

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def get_hyper_pars(X,Y,cov):
        
    gpY = mGP(X,Y,cov)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    res = gpY.min_loglikelihood()
    print(res)
    print(res.x)
    return res.x

#======================================================
datafile = {'OHD51': '/home/ekli/myworks/cosmodata/OHD_51.txt',
        'CC_full': '/home/ekli/myworks/cosmodata/CC_full.dat',
        'OHD_54': '/home/ekli/myworks/cosmodata/OHD_54_CC_BAO.txt',
        'OHD_31': '/home/ekli/myworks/cosmodata/OHD_31_CC.txt',
        'OHD_23': '/home/ekli/myworks/cosmodata/OHD_23_BAO.txt',
        }
H0prior = {'R18': [0.0, 73.52, 1.62],
        'P18': [0.0, 67.27, 0.6],
        }

def add_H0prior(Xs,Ys,covs,which_prior):
    z0, h0, sig0 = H0prior.get(which_prior)
    n = Xs.shape[0]
    X = np.zeros(n+1)
    Y = np.zeros(n+1)
    cov = np.zeros((n+1,n+1))
    X[0] = z0
    X[1:] = Xs
    Y[0] = h0
    Y[1:] = Ys
    cov[0,0] = sig0**2
    cov[1:,1:] = covs
    return X,Y,cov

def get_sigf_lenf(which_dat, which_prior=None):
    dat_file = datafile.get(which_dat)
    Xs,Ys,covs = readin_dat(dat_file)
    if(which_prior is not None):
        X,Y,cov = add_H0prior(Xs,Ys,covs,which_prior)
    else:
        X=Xs
        Y=Ys
        cov=covs
    sigf, lenf = get_hyper_pars(X,Y,cov)
    return X, Y, cov, sigf, lenf

#=======================================================
def rec_hz(zz,X,Y,cov,sigf,lenf):
    gp = fs8GP(X,Y,cov,zz,sigf,lenf)
    hz = gp.rec_mu_arr()
    #cov = gp.cov_over_mumu_mat()
    hcov = gp.cov_mu_mat()
    Ez = gp.rec_mu_over_mu0_arr()
    Ecov = gp.cov_mu_over_mu0_mat()
    return hz, hcov, Ez, Ecov

def rec_fs8(zz,X,Y,cov,sigf,lenf, sig8=0.8, omm=0.3, gam=6/11.0):
    gp = fs8GP(X,Y,cov,zz,sigf,lenf)
    h0 = gp.rec_mu_x(0.0)
    func = lambda x: (omm*(1+x)**3/(gp.rec_mu_x(x)/h0)**2)**gam
    Dint = lambda x: quad(func, 0, x)[0]
    ff = lambda x: sig8*func(x)*np.exp(-Dint(x))
    fs8 = np.array([ff(zi) for zi in zz])
    
    Dints = lambda x: qgaus(func, 0, x)
    ffs = lambda x: sig8*func(x)*np.exp(-Dint(x))
    fs8s = np.array([ff(zi) for zi in zz])
    return fs8, fs8s
    
def rec_theory_fs8(zz,X,Y,cov,sigf,lenf, sig8=0.8, omm=0.3, gam=6/11.0):
    gp = fs8GP(X,Y,cov,zz,sigf,lenf)
    gp.initial_theory_fs8_noH0(sig8,omm,gam)

    fs8 = gp.return_fs8_theory()
    fs8cov = gp.return_fs8cov_theory()
    return fs8, fs8cov

def _lcdm_hz(z, omm=0.3, h0=70):
    return h0*np.sqrt(omm*(1+z)**3 +1-omm)

def lcdm_fs8(zz, sig8=0.8, gam=6.0/11, omm=0.3, h0=70):
    func = lambda x: (omm*h0**2*(1+x)**3/_lcdm_hz(x, omm, h0)**2 )**gam
    Dint = lambda x: quad(func, 0, x)[0]
    ff = lambda x: sig8*func(x)*np.exp(-Dint(x))
    fs8 = np.array([ff(zi) for zi in zz])
    return fs8

#======================================================

def compare_diff_priors_on_H0():
    zz = np.linspace(0,2.1,51)
    
    #dats = ['OHD51','CC_full']
    dats = ['OHD_31', 'OHD_23', 'OHD_54']
    priors = [None, 'R18', 'P18']

    out_dat = {}

    
    for dat in dats:
        out_dat[dat] = {}
        out_dat[dat]['zz'] = zz
        for prior in priors:
            print("="*70)
            print("Start initial out_dat")
            st = time.time()

            X,Y,cov,sigf,lenf = get_sigf_lenf(dat,prior)
            if prior is None:
                prior = 'None'
            hz, hcov, Ez, Ezcov = rec_hz(zz,X,Y,cov,sigf,lenf)
            sig = np.sqrt(np.diag(hcov))

            fs8, fs8s = rec_fs8(zz,X,Y,cov,sigf,lenf)

            theory_fs8, theory_fs8cov = rec_theory_fs8(zz,X,Y,cov,sigf,lenf)

            ###sig0 = np.sqrt(np.diag(hcov))
            ###sig = np.array([h*s for h,s in zip(hz, sig0)])

            ###Ez = np.array([hi/hz[0] for hi in hz])

            ###Ezcov = np.array([[Ez[i]*Ez[j]*(hcov[i,j] -hcov[i,0] -hcov[0,j] + hcov[0,0]) 
            ###    for i in range(hz.shape[0])] for j in range(hz.shape[0])])
            ###Esig2 = np.array([Ei*Ei *sigi +Ei*Ei*hcov[0,0] for Ei,sigi in zip(Ez, np.diag(hcov)) ])            
            ###Esig = np.sqrt(Esig2)
            Esig = np.sqrt(np.diag(Ezcov))

            out_dat[dat][prior] = {}
            out_dat[dat][prior]['mean'] = hz
            out_dat[dat][prior]['sig'] = sig
            out_dat[dat][prior]['Ez_mu'] = Ez
            out_dat[dat][prior]['Ez_sig'] = Esig

            out_dat[dat][prior]['fs8'] = fs8
            out_dat[dat][prior]['fs8s'] = fs8s

            out_dat[dat][priors]['theory_fs8'] = theory_fs8
            out_dat[dat][priors]['theory_fs8cov'] = theory_fs8cov
            #out_dat[dat]['X'] = X
            #out_dat[dat]['Y'] = Y
            #out_dat[dat]['cov'] = cov
            #out_dat[dat]['sigf'] = sigf
            #out_dat[dat]['lenf'] = lenf
            
            et = time.time()
            print("End initial out_dat")
            print("Cost time: %s"%(et - st))
            print("$"*70)
    
    np.save('out_dat', out_dat)

    return out_dat

#======================================================
def test_rec_hz():
    ndim = 41
    zz = np.linspace(0,2.0,ndim)    
    dats = ['OHD_31', 'OHD_23', 'OHD_54', 'OHD51']
    priors = [None, 'R18', 'P18']

    dat = dats[2]
    prior = priors[0]
    
    X,Y,cov,sigf,lenf = get_sigf_lenf(dat,prior)
    hz, hcov, Ez, Ezcov = rec_hz(zz,X,Y,cov,sigf,lenf)

    N = X.shape[0]
    filename = bytes(datafile.get(dat), "utf8")

    status = initial_gapp(N,filename,sigf,lenf)
    status = setup_gapp()
    
    for i in range(ndim):
        print("%5.2f %9.3f %9.6f "%(zz[i], hz[i], np.sqrt(np.diag(hcov))[i]), end='')
        zi = zz[i]
        print("%9.3f %9.6f "%(rec_mu(zi), np.sqrt(rec_covariance(zi, zi)) ), end='')
        print("")

    hz = np.array([rec_mu(z) for z in zz])
    sig = np.array([np.sqrt(rec_covariance(z,z)) for z in zz])
    
    plt.plot(zz, hz, '--b')
    plt.fill_between(zz, hz+sig, hz-sig, color='blue', alpha=0.5)


    plt.show()

#=====================================================
def test_rec_hz_gapp():
    ndim = 41
    zz = np.linspace(0,2.0,ndim)

    dats = ['OHD_31', 'OHD_23', 'OHD_54', 'OHD51']
    priors = [None, 'R18', 'P18']

    for dat in dats:
        #dat = dats[2]
        prior = priors[0]
        
        X,Y,cov,sigf,lenf = get_sigf_lenf(dat,prior)
        
        N = X.shape[0]
        filename = bytes(datafile.get(dat), "utf8")
        
        status = initial_gapp(N,filename,sigf,lenf)
        status = setup_gapp()
        
        h0 = rec_mu(0.0)
        omm = 0.3
        ommh2 = omm*(h0/100)**2
        sig8 = 0.8
        gam = 6.0/11
        status = initial_lcdm(omm,h0)

        #hz = np.array([rec_mu_over_mu0(z) for z in zz])
        #sig = np.array([np.sqrt(rec_cov_mu_over_mu0(z,z)) for z in zz])
        hz = np.array([rec_fs8_with_Hz(z,sig8,ommh2,gam) for z in zz])
        sig = np.array([np.sqrt(rec_cov_fs8_with_Hz(z,z,sig8,ommh2,gam)) for z in zz])
        plt.figure()

        plt.plot(zz, hz, '--b')
        plt.fill_between(zz, hz+sig, hz-sig, color='blue', alpha=0.5)        

        #lhz = np.array([lcdm_hz(z)/h0 for z in zz])
        lhz = np.array([lcdm_growth(z,sig8,gam) for z in zz])
        plt.plot(zz, lhz, ':k')


    plt.show()

#===================================
if __name__ == "__main__":

    test_rec_hz_gapp()
    
    #dats = ['OHD_31', 'OHD_23', 'OHD_54', 'OHD51']
    #priors = [None, 'R18', 'P18']

    #dat = dats[3]
    #prior = priors[0]
    #
    #X,Y,cov,sigf,lenf = get_sigf_lenf(dat,prior)
    ##test_rec_hz()
    #
    #N = X.shape[0]
    ##filename = bytes("/home/ekli/myworks/cosmodata/OHD_51.txt", "utf8")
    #filename = bytes(datafile.get(dat), "utf8")
    ##sigf=170.75260325
    ##lenf=2.61947963

    #status = initial_gapp(N,filename,sigf,lenf)
    #print(status)

    #status = setup_gapp()
    #print(status)
    #
    #zz = np.linspace(0,2,41)
    #
    #hz = np.array([rec_mu(z) for z in zz])
    #sig = np.array([rec_covariance(z,z) for z in zz])
    #
    #plt.plot(zz, hz, '--b')
    #plt.fill_between(zz, hz+sig, hz-sig, color='blue', alpha=0.5)


    #plt.show()

