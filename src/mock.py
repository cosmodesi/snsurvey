#!/usr/bin/env python

import h5py
import matplotlib.pyplot as plt
import numpy
import scipy.stats

fname = '/project/projectdirs/desi/mocks/bgs/MXXL/desi_footprint/v0.0.4/BGS_r20.0.hdf5'
f = h5py.File(fname, "r")

w=numpy.logical_and(f['Data']['z_cos'].value < 0.004, f['Data']['z_cos'].value > 0.002)
abs_mag = f['Data']['abs_mag'].value
weight = numpy.power(-abs_mag[w]/2.5,10)
cumcount, lowerlimit, binsize, _ = scipy.stats.cumfreq(abs_mag[w],weights=weight,numbins=20)
cumcount = cumcount / cumcount[-1]

cumcount = numpy.insert(cumcount,0,0)
zedge = lowerlimit + numpy.arange(len(cumcount))*binsize

zbins = numpy.arange(0,0.5,0.02)
ans = []
ans2 = []
for i in xrange(len(zbins)-1):
    w = numpy.logical_and(f['Data']['z_cos'].value > zbins[i], f['Data']['z_cos'].value < zbins[i+1])
    w1 = numpy.logical_and(w,f['Data']['app_mag'].value <19.5)
    lim  = abs_mag[w1].max()
    ans.append(numpy.interp(lim, zedge,cumcount))
    w1 = numpy.logical_and(w,f['Data']['app_mag'].value <20.)
    lim  = abs_mag[w1].max()
    ans2.append(numpy.interp(lim, zedge,cumcount))


ans=numpy.array(ans)
ans2=numpy.array(ans2)
diff = ans2-ans

ans = ans*.97*.92
ans2 = ans + diff*.92*.77
# plt.plot(zbins[1:],ans,label=r'$r< 19.5$')
plt.plot(zbins[1:],ans2)
# plt.legend()
plt.xlabel('redshift')
plt.ylabel(r'fraction of SNe Ia in galaxies with BGS redshift')
plt.savefig('bg_frac.pdf')