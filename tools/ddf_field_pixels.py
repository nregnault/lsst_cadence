#!/usr/bin/env python

import os
import os.path as op
import sys
import numpy as np
import pylab as pl
import argparse
import logging
logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', 
                    level=logging.INFO)
from ubercal import radec2pix


def ddf_field_radec(cad, proposal_id=5):
    idx = cad['proposalId'] == proposal_id
    logging.info('%s observations with proposal_id=%d' % (idx.sum(), proposal_id))
    cad = cad[idx]
    ra, dec = cad['fieldRA'], cad['fieldDec']
    radec = np.rec.fromarrays((ra,dec), names=['ra','dec'])
    radec = np.unique(radec)
    logging.info('%d DDF fields identified.' % len(radec))
    return radec

def ddf_field_pixels(radec, nside):
    pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='dump central pixel number of DDF fields')
    parser.add_argument('cadence', type=str,
                        help='cadence file')
    parser.add_argument('-n', '--npixels-per-field', 
                        help='number of pixels to dump per field')
    parser.add_argument('--nside', default=64, type=int,
                        help='healpix nside')
    parser.add_argument('--ddf-proposal-id', default=5, 
                        help='DDF proposal ID (not the same for all cadences)')
    args = parser.parse_args()
    print args 
    
    logging.info('loading cadence: %s' % args.cadence)
    cadence = np.load(args.cadence)
    radec = ddf_field_radec(cadence, proposal_id=args.ddf_proposal_id)
    
    r = [radec2pix(ra,dec,nside=args.nside) for (ra,dec) in radec]
    logging.info('pixels: %r' % r)
    
    
