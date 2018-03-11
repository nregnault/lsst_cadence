"""
"""


def radec2pix(ra_deg, dec_deg, nside=1024):
    phi   = ra_deg * RADIAN_PER_DEGREE
    theta = 0.5*np.pi - dec_deg * RADIAN_PER_DEGREE
    return hp.ang2pix(nside, theta, phi, nest=True)


def get_flux_standards(nside, random=0):
    if random <= 0.:
        ra = np.array([0., 12., 32., 55.])
        dec = np.array([0., -22., -12., -45.])
    else:
        ra = np.random.uniform(0., 350., random)
        dec = np.random.uniform(-55, -5., random)
    return radec2pix(ra, dec, nside=nside)
