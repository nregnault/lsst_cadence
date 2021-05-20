"""
"""

import numpy as np
import healpy as hp

RADIAN_PER_DEGREE = np.pi / 180.
DEGREE_PER_RADIAN = 180. / np.pi


def radec2pix(ra_deg, dec_deg, nside=1024):
    """just convert ra,dec into radians and call hp.ang2pix    

    args:
      ra_deg (float): right ascension (degrees)
      dec_deg (float): declination (degrees)
      nside (int): healpix nside

    returns:
      (int) pixel number 
    """
    phi   = ra_deg * RADIAN_PER_DEGREE
    theta = 0.5*np.pi - dec_deg * RADIAN_PER_DEGREE
    return hp.ang2pix(nside, theta, phi, nest=True)



DDF_coordinates = np.rec.fromrecords([('COSMOS',   150.36,    2.840),
                                      ('XMM-LSS',   34.39,   -5.090),   
                                      ('CDFS',      53.00,   -27.44),  
                                      ('ELAIS-S1',   0.,     -45.52),  
                                      ('SPT_DEEP', 349.39,   -63.32),  
                                      ('DDF_820',  119.55,   -43.37),  
                                      ('DDF_858',  187.62,   -42.49),  
                                      ('DDF_1200', 176.63,   -33.15),  
                                      ('DDF_2689', 201.85,     0.93),], names=['name', 'ra', 'dec'])


def get_flux_standards(nside, equatorial=False, ddf_4=False, ddf_all=False, random=0):
    if equatorial:
        ra = np.array([0., 12., 32., 55.])
        dec = np.array([5., -12., -5., 2.2])
    elif ddf_4:
        ra = DDF_coordinates['ra'][:4]
        dec = DDF_coordinates['dec'][:4]
    elif ddf_all:
        ra = DDF_coordinates['ra']
        dec = DDF_coordinates['dec']        
    elif random>0:
        ra = np.random.uniform(0., 350., random)
        dec = np.random.uniform(-55, -5., random)
    else:
        ra = np.array([0., 12., 32., 55.])
        dec = np.array([0., -22., -12., -45.])

    return radec2pix(ra, dec, nside=nside)

