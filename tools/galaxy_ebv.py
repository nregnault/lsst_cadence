#!/usr/bin/env python

"""
Determines and dump a galaxy extinction mask. 
Uses the Schlegel and Finkbeiner maps and the python module 

 git@github.com:kbarbary/sfdmap.git 

"""

import os
import os.path as op
import sys
import sfdmap
import numpy as np
import pylab as pl
import healpy as hp


def main(nside=64, plot=False, ebv_cut=0.2):
    npix = hp.nside2npix(nside)
    
    # pixel coordinates
    ipix = np.arange(npix)
    theta,phi = hp.pixelfunc.pix2ang(nside,ipix,nest=1)
    ra, dec = np.degrees(2*np.pi-phi), -np.degrees(theta-np.pi/2.)

    # E(B-V) values
    m = np.zeros(npix)    
    m[:] = sfdmap.ebv(ra, dec, unit='degree')
    mask = m.copy()

    if ebv_cut>0.:
        m[m>ebv_cut] = hp.UNSEEN
        mask[mask>ebv_cut] = 0.
        mask[mask!=0.] = 1.
    if plot:
        hp.mollview(m, nest=1, min=0., max=0.25)
        hp.graticule()
        hp.mollview(mask, nest=1, min=0., max=2)
        hp.graticule()

    return mask


if __name__ == "__main__":
    mask = main(nside=64, ebv_cut=0.25, plot=1)
    np.save('ebv_mask.npy', mask)

    
    
    

    
    






# checkout this page: 
# https://stackoverflow.com/questions/29702010/healpy-pix2ang-convert-from-healpix-index-to-ra-dec-or-glong-glat



"""
First of all the method pix2ang(nside,indx) gives you the coordinates of pixel with number indx. The pixel number is not directly related to a coordinate, i.e. two consecutive pixel numbers are not necessarily next to each other.

Second, as written in the manual of Healpix (which is the underlying code for healpy) (http://healpix.sourceforge.net/html/csubnode2.htm) the angle theta is defined in range [0,pi] and therefore it cannot directly represent declination [-pi/2,pi/2].

So what I do is I define a transformation and I implement it in two functions, for instance:

def IndexToDeclRa(index):
    theta,phi=hp.pixelfunc.pix2ang(NSIDE,index)
    return -np.degrees(theta-pi/2.),np.degrees(pi*2.-phi)

def DeclRaToIndex(decl,RA):
    return hp.pixelfunc.ang2pix(NSIDE,np.radians(-decl+90.),np.radians(360.-RA))

then the map itself will not be in Decl&RA but if you stick to using IndexToDeclRa and DeclRaToIndex you'll effecively have what you need.
"""

