import numpy as np
import snsims
import sncosmo
from lsst.sims.catUtils.supernovae import SNObject

# Use snsims to get a set of reasonable SN parameter values
zs = [0.3]
nsne = 1
sp = snsims.SimpleSALTDist(nsne,zs, rng=np.random.RandomState()).paramSamples

# create an snsim object that represents a supernova
sn = SNObject()
sncosmo_params= snsims.SimulationTile.getSNCosmoParamDict(sp.ix[0], sn)
sn.set(**sncosmo_params)

# get the sncosmo model out of the snsim object
model = sn.equivalentSNCosmoModel()
model.set(t0=55000.,mwebv=0.)

# get the magnitude on certain dates
model.bandflux('desi', [54990., 55000., 55020.])