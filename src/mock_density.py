#!/usr/bin/env python

import h5py
import matplotlib.pyplot as plt
import numpy
import scipy.stats
from astropy.cosmology import WMAP9 as cosmo
import scipy.integrate as integrate

fname = '/project/projectdirs/desi/mocks/bgs/MXXL/desi_footprint/v0.0.4/BGS_r20.0.hdf5'
f = h5py.File(fname, "r")

fout=open('BGS_nbar.dat','w')

desisr = 14000.*((numpy.pi/180)**2)

zbins = numpy.arange(0,0.1501,1./300)
for i in xrange(len(zbins)-1):
    # w=numpy.logical_and(f['Data']['z_cos'].value < zbins[i+1], f['Data']['z_cos'].value > zbins[i])
    vol = integrate.quad(lambda z: cosmo.differential_comoving_volume(z).value, zbins[i], zbins[i+1])[0] * desisr / (cosmo.h**3)
    vol = vol/1e6
    f.write("{:8.6f} {}\n".format(zbins[i], w.sum()/vol))

f.close()
fout.close()