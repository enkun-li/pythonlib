#!/home/ekli/.local/anaconda3/bin/python
#-*- coding: utf-8 -*-  
#==================================
# File Name: test_like.py
# Author: ekli
# Mail: ekli_091@mail.dlut.edu.cn  
# Created Time: 2019-10-22 10:22:52
#==================================

import numpy as np
from Model.lcdm_omk import Model_lcdm_omk as lcdm
from Model.lcdm_cgr import Model_lcdm_cgr as lcgr
from Model.cosm_tay import Model_cosm_tay as tayl
from Model.cosm_pade import Model_cosm_pade as pade
from Likelihood.HST import Likelihood_HST as hst
from Likelihood.OHD import Likelihood_OHD as ohd
from Likelihood.JLA import Likelihood_JLA as jla
from Likelihood.SDF import Likelihood_6dF as df6
from Likelihood.MGS import Likelihood_MGS as mgs
from Likelihood.DR12Consensus import Likelihood_DR12Consensus as dr12
from Likelihood.DR11quasars import Likelihood_DR11quasars as dr11
from Likelihood.DR12Lya import Likelihood_DR12Lya as dr12a
from Likelihood.DR14quasars import Likelihood_DR14quasars as dr14
from Likelihood.DR11Lya_auto import Likelihood_DR11_Lya_auto as dr11lya_a
from Likelihood.DR11Lya_xcro import Likelihood_DR11_Lya_xcro as dr11lya_x

def test_hst():
    like = hst('/home/ekli/myworks/cosmodata/HST/HST_Riess2018_Gaia.txt')
    for i in range(10):
        Hth = np.random.normal(70,3)
        print("H = %8.3f  chisq = %12.4e"%(Hth, like.loglikelihood(Hth)))
    return

def test_ohd():
    like = ohd('/home/ekli/myworks/cosmodata/OHD/CC_31.txt')
    for i in range(10):
        Hth = np.random.normal(70,3)
        omm = np.random.normal(0.3,0.1)
        omk = np.random.normal(0.0,0.1)
        cosmo = lcdm(Hth, omm, omk, feedback=0)
        print("H = %8.3f omm = %8.3f omk = %8.3f"
                %(Hth, omm, omk), end=' ')
        print("rs = %8.3f"%(cosmo.fit_rs_drag()), end=' ')
        print("chisq = %12.4e"%like.loglikelihood(cosmo.Hofz))
    return

def test_jla():
    datfile = '/home/ekli/myworks/cosmodata/Pantheon/pantheon_zhel.txt'
    covfile = '/home/ekli/myworks/cosmodata/Pantheon/pantheon_cov.txt'
    like = jla(datfile, covfile)
    for i in range(10):
        Hth = np.random.normal(70,3)
        omm = np.random.normal(0.3,0.1)
        omk = np.random.normal(0.0,0.1)
        cosmo = lcdm(Hth, omm, omk)        
        print("H = %8.3f omm = %8.3f omk = %8.3f"
                %(Hth, omm, omk), end=' ')
        print("chisq = %12.4e"%like.loglikelihood(cosmo.distance_DA))
    return

def test_6df():
    like = df6('/home/ekli/myworks/cosmodata/BAO/bao_6dF.txt')
    for i in range(10):
        Hth = np.random.normal(70,3)
        omm = np.random.normal(0.3,0.1)
        omk = np.random.normal(0.0,0.1)
        cosmo = lcdm(Hth, omm, omk)        
        print("H = %8.3f omm = %8.3f omk = %8.3f"
                %(Hth, omm, omk), end=' ')
        print("chisq = %12.4e"%like.loglikelihood(cosmo.distance_DV))
    return

def test_mgs():
    datfile = '/home/ekli/myworks/cosmodata/BAO/bao_MGS.txt'
    probfile = '/home/ekli/myworks/cosmodata/BAO/sdss_MGS_prob.txt'
    like = mgs(datfile, probfile)
    for i in range(10):
        Hth = np.random.normal(70,3)
        omm = np.random.normal(0.3,0.1)
        omk = np.random.normal(0.0,0.1)
        cosmo = lcdm(Hth, omm, omk)        
        print("H = %8.3f omm = %8.3f omk = %8.3f"
                %(Hth, omm, omk), end=' ')
        print("chisq = %12.4e"%like.loglikelihood(cosmo.distance_DV), end=' ')
        print("chisq = %12.4e"%like.loglikelihood_cal(cosmo.distance_DV))
    return

def test_dr12():
    datfile = '/home/ekli/myworks/cosmodata/BAO/bao_DR12Consensus.txt'
    probfile = '/home/ekli/myworks/cosmodata/BAO/BAO_consensus_covtot_dM_Hz.txt'
    like = dr12(datfile, probfile)
    for i in range(10):
        Hth = np.random.normal(70,3)
        omm = np.random.normal(0.3,0.1)
        omk = np.random.normal(0.0,0.1)
        cosmo = lcdm(Hth, omm, omk)        
        print("H = %8.3f omm = %8.3f omk = %8.3f"
                %(Hth, omm, omk), end=' ')
        print("chisq = %12.4e"%like.loglikelihood(cosmo.Hofz, cosmo.distance_DM))
    return

def test_dr11():
    like = dr11('/home/ekli/myworks/cosmodata/BAO/bao_DR11quasars.txt')
    for i in range(10):
        Hth = np.random.normal(70,3)
        omm = np.random.normal(0.3,0.1)
        omk = np.random.normal(0.0,0.1)
        cosmo = lcdm(Hth, omm, omk)        
        print("H = %8.3f omm = %8.3f omk = %8.3f"
                %(Hth, omm, omk), end=' ')
        print("Hz = %8.3f"%cosmo.Hofz(2.36), end=' ')
        print("DA = %8.3f"%cosmo.distance_DA(2.36), end=' ')
        print("chisq = %12.4e"%like.loglikelihood(cosmo.Hofz, cosmo.distance_DA))
    return

def test_dr11lyaa():
    like = dr11lya_a()
    for i in range(10):
        Hth = np.random.normal(70,3)
        omm = np.random.normal(0.3,0.1)
        omk = np.random.normal(0.0,0.1)
        rs = np.random.normal(148.8,4)
        cosmo = lcdm(Hth, omm, omk)        
        print("H = %8.3f omm = %8.3f omk = %8.3f rs = %8.3f"
                %(Hth, omm, omk, rs), end=' ')
        print("Hz = %8.3f"%(3e5/cosmo.Hofz(2.36)/rs/8.708), end=' ')
        print("DA = %8.3f"%(cosmo.distance_DA(2.36)/rs/11.59), end=' ')
        print("chisq = %12.4e"%like.loglikelihood(cosmo.Hofz, cosmo.distance_DA, rs))
    return

def test_dr11lyax():
    like = dr11lya_x()
    for i in range(10):
        Hth = np.random.normal(70,3)
        omm = np.random.normal(0.3,0.1)
        omk = np.random.normal(0.0,0.1)
        cosmo = lcdm(Hth, omm, omk)        
        print("H = %8.3f omm = %8.3f omk = %8.3f"
                %(Hth, omm, omk), end=' ')
        print("Hz = %8.3f"%cosmo.Hofz(2.36), end=' ')
        print("DA = %8.3f"%cosmo.distance_DA(2.36), end=' ')
        print("chisq = %12.4e"%like.loglikelihood(cosmo.Hofz, cosmo.distance_DA))
    return

def test_dr12a():
    like = dr12a('/home/ekli/myworks/cosmodata/BAO/bao_DR12Lya.txt')
    for i in range(10):
        Hth = np.random.normal(70,3)
        omm = np.random.normal(0.3,0.1)
        omk = np.random.normal(0.0,0.1)
        cosmo = lcdm(Hth, omm, omk)        
        print("H = %8.3f omm = %8.3f omk = %8.3f"
                %(Hth, omm, omk), end=' ')
        print("chisq = %12.4e"%like.loglikelihood(cosmo.Hofz, cosmo.distance_DM))
    return

def test_dr14():
    like = dr14('/home/ekli/myworks/cosmodata/BAO/bao_DR14quasars.txt')
    for i in range(10):
        Hth = np.random.normal(70,3)
        omm = np.random.normal(0.3,0.1)
        omk = np.random.normal(0.0,0.1)
        cosmo = lcdm(Hth, omm, omk)        
        print("H = %8.3f omm = %8.3f omk = %8.3f"
                %(Hth, omm, omk), end=' ')
        print("chisq = %12.4e"%like.loglikelihood(cosmo.distance_DV))
    return

def test_dr12a_lcdm_cgr():
    like = dr12a('/home/ekli/myworks/cosmodata/BAO/bao_DR12Lya.txt')
    for i in range(10):
        H0 = np.random.normal(70,3)
        q0 = np.random.normal(-0.55,0.2)
        j0 = np.random.normal(1,0.2)
        s0 = np.random.normal(-0.35, 1)
        l0 = np.random.normal(3.115, 2)
        omk = np.random.normal(0.0, 0.2)
        cosmo = lcgr(H0,q0,j0,s0,l0,omk)
        print("H=%8.3f q0=%8.3f j0=%8.3f s0=%8.3f l0=%8.3f omk=%8.3f"
                %(H0,q0,j0,s0,l0,omk), end=' ')
        print("chisq = %12.4e"%like.loglikelihood(cosmo.Hofz, cosmo.distance_DM))
    return

def test_dr12a_tayl():
    like = dr12a('/home/ekli/myworks/cosmodata/BAO/bao_DR12Lya.txt')
    for i in range(10):
        H0 = np.random.normal(70,3)
        q0 = np.random.normal(-0.55,0.2)
        j0 = np.random.normal(1,0.2)
        s0 = np.random.normal(-0.35, 1)
        l0 = np.random.normal(3.115, 2)
        omk = np.random.normal(0.0, 0.2)
        cosmo = tayl(H0,q0,j0,s0,l0,omk)
        print("H=%8.3f q0=%8.3f j0=%8.3f s0=%8.3f l0=%8.3f omk=%8.3f"
                %(H0,q0,j0,s0,l0,omk), end=' ')
        print("chisq = %12.4e"%like.loglikelihood(cosmo.Hofz, cosmo.distance_DM))
    return

def test_dr12a_pade():
    like  = dr12a('/home/ekli/myworks/cosmodata/BAO/bao_DR12Lya.txt')
    for i in range(10):
        H0 = np.random.normal(70,3)
        q0 = np.random.normal(-0.55,0.2)
        j0 = np.random.normal(1,0.2)
        s0 = np.random.normal(-0.35, 1)
        l0 = np.random.normal(3.115, 2)
        omk = np.random.normal(0.0, 0.2)
        cosmo = pade(H0,q0,j0,s0,l0,omk)
        print("H=%8.3f q0=%8.3f j0=%8.3f s0=%8.3f l0=%8.3f omk=%8.3f"
                %(H0,q0,j0,s0,l0,omk), end=' ')
        print("chisq = %12.4e"%like.loglikelihood(cosmo.Hofz, cosmo.distance_DM))
    return

#=========================
if __name__ == '__main__':
    #test_hst()
    test_ohd()
    #test_jla()
    #test_6df()
    #test_mgs()
    #test_dr12()
    #test_dr11()
    #test_dr12a()
    #test_dr14()
    #print("="*50)
    #test_dr12a_lcdm_cgr()
    #print("="*50)
    #test_dr12a_tayl()
    #print("="*50)
    #test_dr12a_pade()
    print("="*50)
    #test_dr11lyaa()
    print("="*50)
    #test_dr11lyax()
