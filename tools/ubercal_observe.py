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
# import saunerie.selinv as selinv
# import saunerie.linearmodels as lm

from ubercal import lsst
from ubercal import get_flux_standards


RADIAN_PER_DEGREE = np.pi / 180.
DEGREE_PER_RADIAN = 180. / np.pi


OBS_DTYPE=np.dtype([('expnum', 'i4'), ('cell', 'i2'), 
                    # ('ichip', 'i4'), 
                    # ('iraft', 'i4'), 
                    ('pixel', 'i4'), ('mjd', 'f8'),
                    ('band', '|S3'),
                    ('bandname', '|S15'),
                    ('refstar', 'i2'),
                    ('refcell', 'i2'),
                    ('gridobs', 'i2'),
                    ('sky', 'f8'),
                    ('seeing', 'f8'),
                    ('m5', 'f8'),
                    ('exptime', 'f4'),
                    ('proposal', 'i2')])


def main_observe(log, nside, nx=1, ny=1, fast=False, keys={'proposal': 'proposalId'} ):
    """For each pointing, generate the 

    Args:
     - log (ndarray of floats):
         snsim log. List of all pointings
     - nside (int):
         healxpix  of the sphere. 
    
    Returns: numpy.recarray
         observations of individual healpix pixels
    
    .. note: this function expects Ra,Dec in degrees.  The convention
             has changed recently. Double check your inputs.
    
    .. note: we do not set the reference pixels here anymore. 
    """
    if fast:
        f = lsst.FocalPlaneSimple()
        p = f.pixellize()
    else:
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
        
    n = len(log)
    coords = np.rec.fromarrays((np.arange(n), log['Ra'], log['Dec'], 
                                log['mjd'], log['band'], iexp, gridobs, rotation, 
                                log['sky'], log['seeingFwhmEff'], log['fiveSigmaDepth'], 
                                log['exptime'], log[keys['proposal']]))
    
    for (i,r,d,mjd,band,iexp,gridobs,angle, sky, seeing, m5, exptime, proposal) in coords:
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
            z['band'] = band
            z['sky'] = sky
            z['seeing'] = seeing
            z['m5'] = m5
            z['exptime'] = exptime
            z['proposal'] = proposal
            l.append(z)            
            
    # reference pixels
    # not set here anymore. 
    # reference pixels should be decided
    # at a later stage, in the ubercal analysis
    #    if refpixs is not None:
    #        N = len(refpixs)
    #        z = np.zeros(N, dtype=OBS_DTYPE)
    #        z['pixel'] = refpixs
    #        z['refstar'] = 1
    #        l.append(z)
    
    return m, np.hstack(l)


def select_obs_in_mjd_range(log, band, mjd_min, mjd_max):
    """select the observations in the specified mjd range
    """
    logging.info('selecting observations in band=%s and %f<mjd<%f' % (band, mjd_min, mjd_max))
    if band is not None:
        idx  = log['band'] == band
    else:
        idx  = np.ones(len(log)).astype(bool)
    idx &= (log['mjd'] >= mjd_min)
    idx &= (log['mjd'] <= mjd_max)
    logging.info('done. %d observations selected [%d in log]' % (idx.sum(), len(log)))
    return log[idx]


def year_limits(log, year):
    """compute the mjd_min and mjd_max corresponding to a given survey year
    """
    dt = 365.25
    mjd = np.floor(log['mjd'])
    mjd_min, mjd_max = mjd.min(), mjd.max()
    nb_years = np.ceil((mjd_max-mjd_min)/dt)
    assert(year < nb_years)
    mjd_lims = np.arange(mjd.min(), mjd.max()+dt, dt)
    print mjd_min, mjd_max, mjd_lims
    mjd_min, mjd_max = mjd_lims[year], mjd_lims[year+1]
    return mjd_min, mjd_max


def select_survey_year(log, band, year):
    """select the observations taken during the specified survey year
    """
    year = int(year)
    mjd_min, mjd_max = year_limits(log, year)
    
    logging.info('selecting observations in band=%s for survey year: %d' % (band, year))
    logging.info('                       mjd range=[%f,%f]' % (mjd_min, mjd_max))
    if band is not None:
        idx  = log['band'] == band
    else:
        idx  = np.ones(len(log)).astype(bool)    
    idx &= (log['mjd'] >= mjd_min)
    idx &= (log['mjd'] <  mjd_max)
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
    parser.add_argument('--year', default=None, type=int,
                        help='year of survey')
    parser.add_argument('--band',
                        dest='band', default=None, type=str,
                        help='LSST passband [ugrizy]')
    parser.add_argument('--nside', dest='nside', default=1024, type=int,
                        help='HEALPIX pixellization')
    parser.add_argument('--ncellsx', default=1, 
                        type=int, dest='nx',
                        help='number of cells per ccd in the x-direction')
    parser.add_argument('--ncellsy', default=1, 
                        type=int, dest='ny',
                        help='number of cells per ccd in the y-direction')
    parser.add_argument('--plots', default=None,
                        dest='main_plot',
                        help='prepare and dump control plots')
    parser.add_argument('--fast', default=None,
                        dest='fast',
                        help='use a simplistic focal plane')
    parser.add_argument('log', type=str, default=None,
                        help='observation log')

    
    args = parser.parse_args()
    
    if args.log is None:
        parser.error('no observation log specified.')

    # loading data
    if args.year is None:
        log = select_obs_in_mjd_range(np.load(args.log),
                                      band=args.band,
                                      mjd_min=args.mjd_min, mjd_max=args.mjd_max)
    else:
        log = select_survey_year(np.load(args.log), args.band, year=int(args.year))
    
    keys = {'proposal': 'proposalId'}
    if 'survey_id' in log.dtype.names:
        keys['proposal'] = 'survey_id'
    
    # build the observation log
    logging.info('loop on observations ...')
    m,l = main_observe(log, args.nside,
                       nx=args.nx, ny=args.ny, 
                       fast=args.fast,
                       keys=keys)
    logging.info('loop on observations ... done.')    
    
    logging.info(' -> %s' % args.output)
    np.savez(args.output, m=m, l=l)
    
    if args.main_plot is not None:
        hp.mollview(m, nest=1)
        pl.gcf().savefig(args.main_plot, bbox_inches='tight')
        



#    coords = np.rec.fromarrays((np.arange(n), log['Ra'], log['Dec'], log['mjd'], log['band'], iexp, gridobs, rotation, log['sky'], log['seeing'], log['m5sigmadepth'], log['exptime']))
#    for i,(r,d,mjd,band,iexp,gridobs,angle) in enumerate(coords):

#                       refpixs=refpixs,
#                       refpixs=get_flux_standards(args.nside),

#    parser.add_argument('--norefpixs', default=False,
#                        dest='norefpixs', action='store_true',
#                        help='no reference pixel observations')
#    parser.add_argument('--nfluxstd', default=0, type=int,
#                        dest='nb_flux_standards', 
#                        help='number of flux standards to generate')
#    parser.add_argument('--equatorial_standards', default=False, 
#                        action='store_true',
#                        help='use a bunch of equatorial standards (4 in fact) ')

# if args.norefpixs:
#     refpixs = None
# else:
#     refpixs = get_flux_standards(args.nside, random=args.nb_flux_standards, equatorial=args.equatorial_standards)
#     logging.debug('refpixs: %r' % refpixs)
# logging.info('refpixs: %r' % refpixs)

