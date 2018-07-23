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


def get_flux_standards(nside, equatorial=False, random=0):
    if random <= 0.:
        if equatorial:
            ra = np.array([0., 12., 32., 55.])
            dec = np.array([5., -12., -5., 2.2])
        else:
            ra = np.array([0., 12., 32., 55.])
            dec = np.array([0., -22., -12., -45.])
    else:
        ra = np.random.uniform(0., 350., random)
        dec = np.random.uniform(-55, -5., random)
    return radec2pix(ra, dec, nside=nside)

