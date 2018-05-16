#!/usr/bin/env python


"""
Global control plots for the LSST cadences. 

(Draft, final version will be in module lsst_gaia_ubercal)
"""

import os
import os.path as op
import logging
logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',
                    level=logging.DEBUG)
from exceptions import ValueError
from argparse import ArgumentParser

import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab as pl
pl.interactive(0)
import healpy as hp

from saunerie import psf



def _savefig(fig, filename):
    dirname = op.dirname(filename)
    if not op.isdir(dirname):
        os.makedirs(dirname)
    fig.savefig(filename, bbox_inches='tight')
    

def cadence_stats(cad, bands='ugrizy', title=''):
    r = {}
    print '-' * 50
    print '    %s   ' % title
    print '-' * 50
    for band in bands:
        d = cad[cad['band'] == band]
        open_shutter_time = d['exptime'].sum()
        print '%10s: %7d %7.1f' % (band, len(d), open_shutter_time/3600.)
        r[band] = (len(d), open_shutter_time)
    open_shutter_time = cad['exptime'].sum()
    print '%10s: %7d %7.1f' % ('total', len(cad), open_shutter_time/3600.)
    r['all'] = len(cad), open_shutter_time
    return r


def print_all_cadence_stats(bands='ugrizy', plot_band='all'):
    n = len(CADENCE_SHORT_NAMES)
    x = np.arange(n)
    
    nv,exptot = [], []
    for name,fn in zip(CADENCE_SHORT_NAMES,CADENCE_FILES):
        cad = np.load(fn)
        r = cadence_stats(cad, bands=bands, title=fn)
        print r
        s = r.get(plot_band, (0., 0.))
        nv.append(s[0])
        exptot.append(s[1])

    nv = np.array(nv)
    exptot = np.array(exptot)
        
    if plot_band != '':
        pl.figure()
        ax = pl.subplot(211)
        pl.bar(x,nv/1000.)
        pl.ylabel('# visits (x1000)')
    
        ax = pl.subplot(212, sharex=ax)
        pl.bar(x,exptot/3600.)
        pl.subplots_adjust(hspace=0.01)
        pl.ylabel('total exptime [hours]')
        ax.get_xaxis().set_ticks(x)
        ax.get_xaxis().set_ticklabels(CADENCE_SHORT_NAMES, rotation=45)    

    
def number_of_visits(l, band=None, nside=64):
    """
    Given an observation log, compute the number of visits and return
    it as a healpix map.
    """
    npix = hp.nside2npix(nside)
    ipix = l['pixel'][l['band'] == band] if band is not None else l['pixel']
    m = np.bincount(ipix, minlength=npix).astype(float)
    m[m<=0] = hp.UNSEEN
    return m


def plot_number_of_visits(obslog, nside=64, max_nv=400, output_dir=None):
    maps = {}
    for band in "ugrizy":
        logging.info('processing band: %s' % band)
        maps[band] = number_of_visits(obslog, band=band, nside=nside)
    maps['all'] = number_of_visits(obslog, band=None, nside=nside)

    # plot as maps
    fig = pl.figure(num=1, figsize=(8,6))
    fig.clear()
    fig.suptitle('number of visits (full survey)')
    for i,band in zip([1,2,3,4], "griz"):
        nv = maps[band]
        hp.mollview(nv, min=0, max=max_nv, nest=1, fig=1, sub=(2,2,i), 
                    title='%s [%6.0f]' % (band, np.median(nv[nv>0])))
    _savefig(fig, output_dir + os.sep + 'nvisits.png')
    
    # plot as histograms
    fig = pl.figure(num=2, figsize=(8,6))
    fig.clear()    
    fig.suptitle('number of visits (full survey)')
    pl.suptitle('Number of visits (full survey)')
    for i,band in zip([1,2,3,4], "griz"):
        ax = fig.add_subplot(2,2,i)
        nv = maps[band]
        ax.hist(nv[nv>0], range=(0.,max_nv), bins=100)
        ax.set_title(band)
    _savefig(fig, output_dir + os.sep + 'nvisits_hists.png')
    
    return maps
    

def average_seeing(l, band, nside=64):
    npix = hp.nside2npix(nside)
    idx = l['band'] == band
    ipix = l['pixel'][idx]
    seeing = l['seeing'][idx]
    nv = np.bincount(ipix, minlength=npix)
    seeing = np.bincount(ipix, minlength=npix, weights=seeing)
    seeing = seeing / nv
    seeing[np.isnan(seeing)] = hp.UNSEEN
    return seeing


def plot_average_seeing(obslog, nside=64, output_dir=None):
    fig = pl.figure(num=1, figsize=(8,6))
    fig.clear()
    pl.suptitle('Average seeing (full survey)')
    for i,band in enumerate("griz"):
        seeing = average_seeing(obslog, band, nside=nside)
        hp.mollview(seeing, min=0.5, max=1.2, nest=1, fig=1, sub=(2,2,i+1), title='%s [%5.2f\'\']' % (band, np.median(seeing[seeing>0])))
    _savefig(fig, output_dir + os.sep + 'seeing.png')


def survey_depth(l, band, etc, nside=64):
    npix = hp.nside2npix(nside)
    out = np.zeros(npix)
    
    idx = l['band'] == band
    ipix = l['pixel'][idx]
    m5 = l['m5'][idx]
    
    f5 = etc.mag_to_flux(m5, [etc.instrument.name + '::' + band])
    #                         np.core.defchararray.add(etc.instrument.name + '::', l['band'][idx]))
    var = (f5 / 5.)**2
    stack_variance = np.bincount(ipix, minlength=npix, weights=1./var)
    #    nv = np.bincount(ipix, minlenwgth=npix)
    idx = stack_variance > 0
    out[idx] = 1. / stack_variance[idx]
    out[idx] = 5 * np.sqrt(out[idx])
    out[idx] = etc.flux_to_mag(out[idx], [etc.instrument.name + '::' + band])
    out[~idx] = hp.UNSEEN

    return out


def plot_survey_depth(l, etc, nside=64, maxv=29., output_dir=None):
    maps = {}
    for band in "ugrizy":
        logging.info('processing band: %s' % band)
        maps[band] = survey_depth(l, band, etc, nside=nside)

    # plot as maps
    fig = pl.figure(num=1, figsize=(8,6))
    fig.clear()
    pl.suptitle('survey depth (full survey)')
    for i,band in zip([1,2,3,4], "griz"):
        m = maps[band]
        hp.mollview(m, min=24.5, max=maxv, nest=1, fig=1, sub=(2,2,i),
                    title='%s [%5.2f]' % (band, np.median(m[m>0.])))
    _savefig(fig, output_dir + os.sep + 'survey_depth.png')
        
    # plot as histograms
    fig = pl.figure(num=2, figsize=(8,6))
    fig.clear()
    pl.suptitle('Survey depth (full survey)')
    for i,band in zip([1,2,3,4], "griz"):
        ax = fig.add_subplot(2,2,i)
        nv = maps[band]
        ax.hist(nv[nv>0], range=(24.5,maxv), bins=100, log=True)
        ax.set_title(band)
    _savefig(fig, output_dir + os.sep + 'survey_depth_hist.png')
    
    return maps
    

def average_cadence(l, band=None, nside=64, exclude_bands=[], season_gap=100):
    """
    """
    npix = hp.nside2npix(nside)
    cad = np.zeros(npix)
    
    N = len(l)
    idx_excl = np.zeros(N).astype(bool)
    for bn in exclude_bands:
        idx_excl &= (l['band'] == bn)
    idx_band = l['band'] == band    
    ll = l[~idx_excl & idx_band]

    ipix = np.unique(ll['pixel'])
    mjd = np.floor(ll['mjd'])
    for j,i in enumerate(ipix):
        if j%100 == 0:
            logging.info('pixel: %d / %d' % (j, len(ipix)))
        idx = ll['pixel'] == i
        m = np.unique(mjd[idx])
        dt = m[1:] - m[:-1]
        dt = dt[dt>0.]
        cad[i] = np.mean(dt[dt<season_gap])
        
    cad[cad<=0.] = hp.UNSEEN
    
    return cad
    

def plot_average_cadence(data, nside=64, output_dir=None):
    fig = pl.figure(num=1, figsize=(8,6))
    fig.clear()
    pl.suptitle('median $\Delta t$ between obs (Full survey)')
    # main cadence
    for i,band in zip([1,2,3,4],'griz'):
        logging.info('processing band: %s' % band)
        cad = average_cadence(data, band=band, nside=64)
        hp.mollview(cad, min=0, max=50, nest=1, fig=1, sub=(2,2,i), 
                    title='%s [%5.1f days]' % (band, np.median(cad[cad>0.])))
    _savefig(fig, output_dir + os.sep + 'average_cadence.png')


if __name__ == '__main__':
    parser = ArgumentParser(
        description="""global stats for the cadence specified in argumement""")
    parser.add_argument('-O', '--output_dir',
                        default=None,
                        help='output directory')
    # parser.add_argument('cadence', type=str, default=None,
    #                     help='cadence_filename')
    parser.add_argument('obslog', type=str,
                        help='obslog')
    args = parser.parse_args()
    
    # open the obslog 
    #    cad = np.load(args.cadence_filename)
    with np.load(args.obslog) as f:
        obslog = f['l']
    
    
    # global stats
    #    r = cadence_stats(cad, bands='ugrizy', title=args.cadence_filename)

    # numbers of visits
    plot_number_of_visits(obslog, nside=64, max_nv=400, output_dir=args.output_dir)

    # average seeing
    plot_average_seeing(obslog, nside=64, output_dir=args.output_dir)
    
    # survey depth
    etc = psf.find('LSSTPG')
    plot_survey_depth(obslog, etc, nside=64, maxv=28., output_dir=args.output_dir)

    # ... aaaaand finally ! (long !)
    if True:
        plot_average_cadence(obslog, nside=64, output_dir=args.output_dir)
