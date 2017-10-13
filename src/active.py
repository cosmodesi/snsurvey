#!/usr/bin/env python
import numpy
import scipy.integrate as integrate
import control
import snrates
import matplotlib.pyplot as plt

def integrand(z, rlim):
    # print z, snrates.rodney2014.rate(z), control.controlTime(z,rlim)
    return snrates.rodney2014.rate(z) * \
        control.controlTime(z,rlim) * \
        snrates.cosmo.differential_comoving_volume(z).value/(1+z)/365.25

lmag = numpy.array([19.5,20.5,21.5])
zbinedge = numpy.arange(0,0.3501,.05)

ans=[]
for lm in lmag:
    ans_=[]
    cum=0.
    for i in xrange(1,len(zbinedge)):
        result = snrates.sr_deg * integrate.quad(integrand, zbinedge[i-1], zbinedge[i], args=lm, epsabs=1e-4)[0]
        print lm, zbinedge[i],result
        cum = cum+result
        ans_.append(cum)
    ans_ = numpy.array(ans_)
    ans.append(ans_)



for lm, act in zip(lmag,ans):
    plt.plot(zbinedge[1:],act, label = '$r_{{lim}} = {}$'.format(str(lm)))

plt.legend(loc=2)
plt.xlabel(r'$z$')
plt.ylabel(r'Cumulative $N_{{SN Ia}}$ per s.d.')
plt.savefig('cumulative.pdf')


# cum=0.
# print '{:>6s} {:>6s} {:>10s} {:>10s}'.format('zmin','zmax','# bin','# cum')
# print '{:>33s}'.format('1/(sd yr)')
# for i in xrange(1,len(zbinedge)):
#     result = snrates.sr_deg * integrate.quad(integrand, zbinedge[i-1], zbinedge[i], args=20.)[0]
#     cum = cum+result
#     print '{:>6.2f} {:>6.2f} {:>10.3f} {:>10.3f}'.format(zbinedge[i-1], zbinedge[i], result,cum)