from astropy.cosmology import FlatLambdaCDM
import fitsio
import numpy as np
import sys

def get_derived_quantities(template_path, coeff, z): 
  """
  Given a path to a kcorrect template file, as well as 
  kcorrect coefficients for a set of galaxies and return derived
  quantities
  
  inputs
  ------
  template_path -- str
    path to kcorrect templates
  coeff -- (N x N_templates) array
    array of kcorrect coefficients
  z -- (N,) array
    redshifts of galaxies
    
  returns
  -------
  sfr -- (N,) array
    Array of star formation rates for galaxies
  met -- (N,) array
    Array of metallicities for galaxies
  smass -- (N,) array
    Array of stellar masses for galaxies
  """
  
  #read in relevant template info
  sfh_tot = fitsio.read(template_path, 12)
  sfh_met = fitsio.read(template_path, 13)
  mass_tot = fitsio.read(template_path, 17)
  
  #get angular diameter distances
  cosmo = FlatLambdaCDM(100, 0.286)
  da    = cosmo.angular_diameter_distance(z).value  
  
  smass = np.dot(mass_tot, coeff.T) * (da * 1e6 / 10) ** 2
  ssfr  = np.dot(coeff, sfh_tot)
  met   = np.dot(coeff, sfh_tot * sfh_met) / ssfr
  sfr   = ssfr * smass[:,np.newaxis]
  
  #get the values at z_galaxy
  met   = met[:,-1]
  sfr   = sfr[:,-1]
  
  return sfr, met, smass
  

if __name__=='__main__':
    kfile = sys.argv[1] #name of file containing 
    filename = sys.argv[2] #name of galaxy catalog file

    zmax = 0.1

    galaxies  = fitsio.FITS(filename)[-1].read(columns=['COEFFS','Z']) # read relevant info from files

    filt= np.array([np.isfinite(g['COEFFS']).all() and g['Z'] < zmax for g in galaxies])
    print galaxies.shape
    galaxies = galaxies[filt]
    print galaxies.shape
    sfr, met, smass = get_derived_quantities(kfile,galaxies['COEFFS'],galaxies['Z'])
    print sfr[0], met[0], smass[0]