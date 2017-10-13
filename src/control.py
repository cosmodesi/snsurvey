#!/usr/bin/env python
import numpy
import sncosmo
import scipy.optimize

import matplotlib.pyplot as plt

model=sncosmo.Model(source='salt2-extended')

def f(t ,rlim):
    # print t, model.bandflux('desr',t, zp = rlim, zpsys='ab')
    return model.bandflux('desr',t, zp = rlim, zpsys='ab')-1.

def controlTime(z,rlim):
    model.set(z=z, t0=55000.)
    model.set_source_peakabsmag(absmag=-19.3,band='bessellb',magsys='ab')
    pre = scipy.optimize.fsolve(f, 55000.-15*(1+z) ,args=(rlim),xtol=1e-8)
    post = scipy.optimize.fsolve(f, 55000.+20*(1+z) ,args=(rlim),xtol=1e-8)
    return max(post[0]-pre[0],0)
    # print scipy.optimize.fsolve(f, 55000.+40,args=(rlim),factor=1.,xtol=1e-8)

def plot():

    lmag = numpy.arange(19.5,21.6,0.5)
    zs = numpy.arange(0.02, 0.2501,0.02)

    ans = []
    for lm in lmag:
        ans_=[]
        for z in zs:
            ans_.append(controlTime(z,lm))
        ans.append(ans_)


    for lm, ct in zip(lmag, ans):
        plt.plot(zs, ct, label = '$r_{{lim}} = {}$'.format(str(lm)))

    plt.xlabel(r'$z$')
    plt.ylabel(r'control time (days)')
    plt.legend()
    plt.show()
