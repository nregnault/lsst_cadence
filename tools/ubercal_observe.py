#!/usr/bin/env python 

"""


"""

import os
import os.path as op
import sys
import argparse

import logging
logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', 
                    level=logging.INFO)

import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab as pl
pl.interactive(0)

import numpy as np
import healpy as hp
from mx.DateTime import DateTimeFrom
import croaks
import saunerie
from saunerie.fitparameters import FitParameters 
import saunerie.selinv as selinv
import saunerie.linearmodels as lm

from ubercal import lsst
from ubercal import get_flux_standards


RADIAN_PER_DEGREE = np.pi / 180.
DEGREE_PER_RADIAN = 180. / np.pi

OBS_DTYPE=np.dtype([('expnum', 'i4'), ('cell', 'i2'), 
                    # ('ichip', 'i4'), 
                    # ('iraft', 'i4'), 
                    ('pixel', 'i4'), ('mjd', 'f8'),
                    ('refstar', 'i2'),
                    ('refcell', 'i2'),
                    ('gridobs', 'i2')])


def main_observe(log, nside, refpixs, nx=1, ny=1):
    """For each pointing, generate the 

    Args:
     - log: 
         snsim log. List of all pointings
     - nside: int
         healxpix  of the sphere. 
     - refpixels: list
         coordinates of the reference pixels 

    Returns: numpy.recarray
         observations of individual healpix pixels

    .. note: this function expects Ra,Dec in degrees.  The convention
             has changed recently. Double check your inputs.
    """
    f = lsst.FocalPlane()
    p = f.pixellize(nx,ny)
    l = []
    m = np.zeros(hp.nside2npix(nside))

    n = len(log)
    iexp = np.arange(n, dtype=int)
    if 'gridobs' in log.dtype.names:
        gridobs = log['gridobs']
    else:
        gridobs = np.zeros(n).astype(bool)
    if 'rotation' in log.dtype.names:
        rotation = log['rotation']
    else:
        rotation = np.zeros(n)
    coords = np.vstack((log['Ra'],
                        log['Dec'],
                        log['mjd'],
                        iexp,
                        gridobs,
                        rotation)).T
    for i,(r,d,mjd,iexp,gridobs,angle) in enumerate(coords):
        if i%100 == 0:
            logging.info('exposure: % 6d [RA=% 6.6f DEC=% +6.5f]' % (i,r,d))
        cells = p.copy()
        if np.abs(angle)>0.:
            cells = f.rotate(cells, theta=angle * 180. / np.pi)
        f.to_uv(cells, (r,d))

        # main observations
        for c in cells:
            poly = np.vstack([c['p0'], c['p1'], c['p2'], c['p3']])
            try:
                ipix = hp.query_polygon(nside, poly, nest=True, inclusive=False)
            except:
                continue
            m[ipix] += 1
            N = len(ipix)
            z = np.zeros(N, dtype=OBS_DTYPE)
            z['expnum'] = iexp
            z['cell'] = c['icell']
            z['pixel'] = ipix
            z['mjd'] = mjd
            z['gridobs'] = gridobs
            l.append(z)
            
        
        # reference pixels
        if refpixs is not None:
            N = len(refpixs)
            z = np.zeros(N, dtype=OBS_DTYPE)
            z['pixel'] = refpixs
            z['refstar'] = 1
            l.append(z)

    return m, np.hstack(l)


def select(log, band, mjd_min, mjd_max):
    logging.info('selecting observations in band=%s and %f<mjd<%f' % (band, mjd_min, mjd_max))    
    idx  = log['band'] == band
    idx &= (log['mjd'] >= mjd_min)
    idx &= (log['mjd'] <= mjd_max)
    logging.info('done. %d observations selected [%d in log]' % (idx.sum(), len(log)))
    return log[idx]
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""For a given observation log, create an 'obs-file'.  Each line of
        the obs-file describes the observation of one healpix in one
        focal plane cell.""")
    parser.add_argument('-o', '--output',
                        default='obs.npy',
                        help='output obs file')
    parser.add_argument('--mjd_min',
                        dest='mjd_min', default=59580., type=float,
                        help='observation period start')
    parser.add_argument('--mjd_max',
                        dest='mjd_max', default=59945., type=float,
                        help='observation period end')
    parser.add_argument('--band',
                        dest='band', default='z', type=str,
                        help='LSST passband [ugrizy]')
    parser.add_argument('--nside', dest='nside', default=1024, type=int,
                        help='HEALPIX pixellization')
    parser.add_argument('--ncellsx', default=1, type=int, dest='nx',
                        help='number of cells per ccd in the x-direction')
    parser.add_argument('--ncellsy', default=1, type=int, dest='ny',
                        help='number of cells per ccd in the y-direction')
    parser.add_argument('--plots', default=None,
                        dest='main_plot',
                        help='prepare and dump control plots')
    parser.add_argument('log', type=str, default=None,
                        help='observation log')

    
    args = parser.parse_args()
    
    if args.log is None:
        parser.error('no observation log specified.')

    # loading data
    log = select(np.load(args.log),
                         band=args.band,
                         mjd_min=args.mjd_min, mjd_max=args.mjd_max)

    # build the observation log
    logging.info('loop on observations ...')
    m,l = main_observe(log, args.nside, refpixs=get_flux_standards(args.nside),
                       nx=args.nx, ny=args.ny)
    logging.info('loop on observations ... done.')    
    
    logging.info(' -> %s' % args.output)
    np.savez(args.output, m=m, l=l)
    
    if args.main_plot is not None:
        hp.mollview(m, nest=1)
        pl.gcf().savefig(args.main_plot, bbox_inches='tight')
        
