#!/usr/bin/env python

import numpy
import scipy.stats

logssfr_lim = -11.2

# rate model
# rate = c_M (M-<M>) + c_sSFR (sSFR - <sSFR>) + <rate>

# 0: SN Ia star-forming
# 1: SN Ia passive
# 2: SN II
# 3: null


coeffs  = []

# sSFR

m_mn = 2.3
ssfr_mn = 70.
r_mn = 0.129
m = numpy.array([0.6,2.7,7.7])
r_m = numpy.array([0.3,0.15,0.076])
ssfr =  numpy.array([18,70,200.])
r_ssfr = numpy.array([0.086,0.129,0.265])
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
c_M, _, _, _, _ = scipy.stats.linregress(m-m_mn, r_m-r_mn)
c_sSFR, _, _, _, _ = scipy.stats.linregress(ssfr-ssfr_mn, r_ssfr-r_mn)

coeffs.append((m_mn, ssfr_mn, r_mn, c_M, c_sSFR))
coeffs.append((0,0,0,0,0.))

coeffs = numpy.array(coeffs, order='F')

# units
# mass : 1e10 Msun
# sSFR : 1e-12 / yr
# rate : 1e-12 / yr/ M_sun

def SNIa(mass, ssfr):
    ind = numpy.log10(numpy.array(ssfr,copy=False)*1e-12) < logssfr_lim    #if true 1 is passive, 0 is ssfr
    ind = ind.astype(int)

    return coeffs[ind,2] + coeffs[ind,3]*(numpy.array(mass,copy=False)-coeffs[ind,0])+ \
        coeffs[ind,4]*(numpy.array(ssfr,copy=False)-coeffs[ind,1]) 

def SNII(mass, ssfr):
    ind = numpy.log10(numpy.array(ssfr,copy=False)*1e-12) < logssfr_lim    #if true 1 is passive, 0 is ssfr
    ind = ind.astype(int)+2

    return coeffs[ind,2] + coeffs[ind,3]*(numpy.array(mass,copy=False)-coeffs[ind,0])+ \
        coeffs[ind,4]*(numpy.array(ssfr,copy=False)-coeffs[ind,1]) 



# print SNIa([1.2,2.2], [200,68])
# print SNIa(1.2, 200)

# print SNII([1.,1.],[170.,1.])
# print SNII(1., 170)