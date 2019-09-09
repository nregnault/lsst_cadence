#!/usr/bin/env python 

import os
import os.path as op
import sys
import numpy as np 
import healpy as hp
import argparse

from ubercal import hpixlog


def process_obslog(obslog, ddf_id=5):
    _, _, ddf_pixels = hpixlog.extract_ddf_data(obslog, ddf_id=ddf_id)
    ddf_pixels.sort()
    print len(ddf_pixels)
    u = ddf_pixels[1:] - ddf_pixels[:-1]
    print u.min(), u.max(), u
    idx = np.abs(u)>1000
    print np.sum(idx), ' clusters identified'
    bb = ddf_pixels[np.argwhere(idx)]
    boundaries = np.vstack((bb[:-1], bb[1:]))
    print boundaries
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="extract N DDF pixels in each of the DDF identified in that cadence")
    parser.add_argument('--ddf_id',
                        default=5, type=int, 
                        help='DDF proposal id for this cadence (5, generally, 1 in genral for altsched)')
    parser.add_argument('--nb_stars_per_ddf',
                        default=1, type=int, 
                        help='how many pixels should I choose from each DDF')
    parser.add_argument('obslog', type=str, default=None, 
                        help='observing log')
    args = parser.parse_args()
    
    print args 
    
    obslog = None
    with np.load(args.obslog) as f:
        obslog = f['l']
    process_obslog(obslog, ddf_id=args.ddf_id)

    


