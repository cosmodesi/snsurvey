#!/usr/bin/env python
import numpy
import scipy.integrate as integrate
import scipy.stats
from astropy.cosmology import WMAP9 as cosmo
import matplotlib.pyplot as plt

sr_deg = (numpy.pi/180)**2

class graur15:

    global logssfr_lim
    global coeffs

    logssfr_lim = -11.2

    # rate model
    # rate = c_M (M-<M>) + c_sSFR (sSFR - <sSFR>) + <rate>

    # 0: SN Ia star-forming
    # 1: SN Ia passive
    # 2: SN II
    # 3: null


    coeffs  = []

    
        # Input units from table

        # Stellar mass    1e10 solar mass
        # sSFR            1e-12 per yr
        # SN rate         1e-12 per yr pr solar mass
    


    # sSFR

    m_mn = 2.3
    ssfr_mn = 70.
    r_mn = 0.129
    m = numpy.array([0.6,2.7,7.7])
    r_m = numpy.array([0.3,0.15,0.076])
    ssfr =  numpy.array([18,70,200.])
    r_ssfr = numpy.array([0.086,0.129,0.265])

    m_mn = numpy.log10(m_mn)+10
    m = numpy.log10(m)+10

    ssfr_mn = numpy.log10(ssfr_mn) -12
    ssfr = numpy.log10(ssfr) -12

    r_mn = numpy.log10(r_mn) - 12
    r_m = numpy.log10(r_m) -12
    r_ssfr = numpy.log10(r_ssfr) -12 

    c_M, _, _, _, _ = scipy.stats.linregress(m-m_mn, r_m-r_mn)
    c_sSFR, _, _, _, _ = scipy.stats.linregress(ssfr-ssfr_mn, r_ssfr-r_mn)

    coeffs.append((m_mn, ssfr_mn, r_mn, c_M, c_sSFR))

    # passive

    m_mn = 10.
    ssfr_mn = 0.8
    r_mn = 0.081
    m = numpy.array([3.6,7.8,18])
    r_m = numpy.array([0.091,0.126,0.056])
    ssfr =  numpy.array([0.39,0.81,2.2])
    r_ssfr = numpy.array([0.068,0.107,0.078])

    m_mn = numpy.log10(m_mn)+10
    m = numpy.log10(m)+10

    ssfr_mn = numpy.log10(ssfr_mn) -12
    ssfr = numpy.log10(ssfr) -12

    r_mn = numpy.log10(r_mn) - 12
    r_m = numpy.log10(r_m) -12
    r_ssfr = numpy.log10(r_ssfr) -12 

    c_M, _, _, _, _ = scipy.stats.linregress(m-m_mn, r_m-r_mn)
    c_sSFR, _, _, _, _ = scipy.stats.linregress(ssfr-ssfr_mn, r_ssfr-r_mn)

    coeffs.append((m_mn, ssfr_mn, r_mn, c_M, c_sSFR))


    # SNe II

    m_mn = 1.
    ssfr_mn = 130.
    r_mn = 0.52
    m = numpy.array([0.08,0.35,2.1])
    r_m = numpy.array([5.5,1.8,0.187])
    ssfr =  numpy.array([69.,170,400])
    r_ssfr = numpy.array([0.25,0.66,2.62])

    m_mn = numpy.log10(m_mn)+10
    m = numpy.log10(m)+10

    ssfr_mn = numpy.log10(ssfr_mn) -12
    ssfr = numpy.log10(ssfr) -12

    r_mn = numpy.log10(r_mn) - 12
    r_m = numpy.log10(r_m) -12
    r_ssfr = numpy.log10(r_ssfr) -12 
    c_M, _, _, _, _ = scipy.stats.linregress(m-m_mn, r_m-r_mn)
    c_sSFR, _, _, _, _ = scipy.stats.linregress(ssfr-ssfr_mn, r_ssfr-r_mn)

    coeffs.append((m_mn, ssfr_mn, r_mn, c_M, c_sSFR))
    coeffs.append((0,0,0,0,0.))

    coeffs = numpy.array(coeffs, order='F')


    # units
    # mass : 1e10 Msun
    # sSFR : 1e-12 / yr
    # rate : 1e-12 / yr/ M_sun
    @staticmethod
    def SNIa(mass, ssfr):
        ind = numpy.log10(numpy.array(ssfr,copy=False)) < logssfr_lim    #if true 1 is passive, 0 is ssfr
        ind = ind.astype(int)

        return 10**(coeffs[ind,2] + coeffs[ind,3]*(numpy.array(numpy.log10(mass),copy=False)-coeffs[ind,0])+ \
            coeffs[ind,4]*(numpy.array(numpy.log10(ssfr),copy=False)-coeffs[ind,1]))
    
    @staticmethod
    def SNII(mass, ssfr):
        ind = numpy.log10(numpy.array(ssfr,copy=False)) < logssfr_lim    #if true 1 is passive, 0 is ssfr
        ind = ind.astype(int)+2
        ans = 10**(coeffs[ind,2] + coeffs[ind,3]*(numpy.array(numpy.log10(mass),copy=False)-coeffs[ind,0])+ \
            coeffs[ind,4]*(numpy.array(numpy.log10(ssfr),copy=False)-coeffs[ind,1]))
        if isinstance(ind, numpy.int64) and ind==3:
            ans = 0
        if not isinstance(ind, numpy.int64):
            ans[ind ==3]=0
        return ans


class dilday2010:

    # http://iopscience.iop.org/article/10.1088/0004-637X/713/2/1026/meta
    # rV = (2.69+0.34+0.21 -0.30-0.01)x10-5 SNe yr-1 Mpc-3 (H 0/(70 km s-1 Mpc-1))3
    global restrate
    restrate =  2.69e-5*(cosmo.h/0.70)**3

    @staticmethod
    def tabulate(zbinedge):
        cum=0.
        for i in xrange(1,len(zbinedge)):
            result = restrate * sr_deg * integrate.quad(lambda z: \
                cosmo.differential_comoving_volume(z).value/(1+z), zbinedge[i-1], zbinedge[i])[0]
            cum = cum+result
            print zbinedge[i-1], zbinedge[i], result,cum

class rodney2014:

    # https://arxiv.org/pdf/1401.7978.pdf
    # SNuVol = 10-4 yr-1 Mpc-3 h703

    global redshift
    global snu

    redshift  = numpy.array([0,0.25,0.75,1.25,1.75,2.25])
    snu  = 1e-4*numpy.array([0.269,0.36,0.51,.64,.72,.49])

    @staticmethod
    def rate(z):
        return ((cosmo.h/.7)**3)*numpy.interp(z,redshift,snu)

    @staticmethod
    def integrand(z):
        return ((cosmo.h/.7)**3)*numpy.interp(z,redshift,snu) *\
            cosmo.differential_comoving_volume(z).value/(1+z)

    @staticmethod
    def tabulate(zbinedge):
        cum=0.
        print '{:>6s} {:>6s} {:>10s} {:>10s}'.format('zmin','zmax','# bin','# cum')
        print '{:>33s}'.format('1/(sd yr)')
        for i in xrange(1,len(zbinedge)):
            result = sr_deg * integrate.quad(rodney2014.integrand, zbinedge[i-1], zbinedge[i])[0]
            cum = cum+result
            print '{:>6.2f} {:>6.2f} {:>10.3f} {:>10.3f}'.format(zbinedge[i-1], zbinedge[i], result,cum)


    @staticmethod
    def plot(zbinedge):
        ans = []
        cum=0.
        for i in xrange(1,len(zbinedge)):
            result = sr_deg * integrate.quad(rodney2014.integrand, zbinedge[i-1], zbinedge[i])[0]
            cum = cum+result
            ans.append(cum)
        ans= numpy.array(ans)
        plt.plot(zbinedge[1:],ans)
        plt.xlabel(r'$z$')
        plt.yscale("log", nonposy='clip')
        plt.ylabel(r'Cumulative $N_{{SN Ia}}$ per s.d. per yr')
        plt.savefig('total.pdf')


# dilday2010.tabulate(numpy.arange(0,.5,.1))
# rodney2014.tabulate(numpy.arange(0,1.2,.05))



# print graur15.SNIa(2.7e10,70e-12)
# print graur15.SNIa(7.8e10,.81e-12)
# print graur15.SNII(0.35e10,170e-12)
# print graur15.SNII(0.35e10,1e-12)