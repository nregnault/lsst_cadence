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
from saunerie.robuststat import mad
from ubercal import hpixlog

CADENCE_FILES = ['alt_sched.npy',  'alt_sched_rolling.npy',  'feature_baseline_10yrs.npy',
                 'feature_rolling_half_mask_10yrs.npy',  'feature_rolling_twoThird_10yrs.npy',
                 'minion_1016.npy']
CADENCE_SHORT_NAMES = ['AltSched',  'AltSchedRolling',  'FeatureBaseline',
                       'FeatureRolling1/2',  'FeatureRolling2/3',
                       'Minion']


def _savefig(fig, filename):
    dirname = op.dirname(filename)
    if not op.isdir(dirname):
        os.makedirs(dirname)
    fig.savefig(filename, bbox_inches='tight')


def cadence_stats(cad, bands='ugrizy', title='', silent=False):
    r = {}

    open_shutter_time = cad['exptime'].sum()
    nvisits_tot = len(cad)
    r['all'] = nvisits_tot, open_shutter_time
    
    header, nvisits, ostime, prct = '', '', '', ''
    for band in bands:
        header += '%10s & ' % band
    header += ' total '
    
    for band in bands:
        d = cad[cad['band'] == band]
        open_shutter_time = d['exptime'].sum()
        nvisits += '%7d  & ' % len(d)
        ostime += '%7.1f & ' % (open_shutter_time / 3600.)
        prct += '%7.1f & ' % (100. * len(d) / nvisits_tot)
        r[band] = (len(d), open_shutter_time)
    
    nvisits += '%d' % len(cad)
    ostime += '%8.1f' % (open_shutter_time / 3600.)
    
    if not silent:
        print '-' * 50
        print '    %s   ' % title
        print '-' * 50
        print header
        print nvisits
        print prct
        print ostime
        
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


def plot_number_of_visits(obslog, nside=64, max_nv=400, output_dir=None, savemaps=False, filename='nvisits.npz', plotname='nvisits.png'):
    maps = {}
    for band in "ugrizy":
        logging.info('processing band: %s' % band)
        maps[band] = number_of_visits(obslog, band=band, nside=nside)
    maps['all'] = number_of_visits(obslog, band=None, nside=nside)
    
    # save healpix maps 
    if savemaps:
        np.savez(output_dir + os.sep + filename, 
                 **maps)

    # plot as maps
    fig = pl.figure(num=1, figsize=(8,6))
    fig.clear()
    fig.suptitle('number of visits (full survey)')
    for i,band in zip([1,2,3,4], "griz"):
        nv = maps[band]
        hp.mollview(nv, min=0, max=max_nv, nest=1, fig=1, sub=(2,2,i), 
                    title='%s [%6.0f]' % (band, np.median(nv[nv>0])))
    _savefig(fig, output_dir + os.sep + plotname)
    
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


def plot_average_seeing(obslog, nside=64, output_dir=None, savemaps=False):
    maps = {}
    fig = pl.figure(num=1, figsize=(8,6))
    fig.clear()
    pl.suptitle('Average seeing (full survey)')
    for i,band in enumerate("griz"):
        seeing = average_seeing(obslog, band, nside=nside)
        hp.mollview(seeing, min=0.5, max=1.2, nest=1, fig=1, sub=(2,2,i+1), title='%s [%5.2f\'\']' % (band, np.median(seeing[seeing>0])))
        maps[band] = seeing
    _savefig(fig, output_dir + os.sep + 'seeing.png')
    if savemaps:
        np.savez(output_dir + os.sep + 'seeing.npz', 
                 **maps)


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


def plot_survey_depth(l, etc, nside=64, maxv=29., output_dir=None, savemaps=False):
    maps = {}
    for band in "ugrizy":
        logging.info('processing band: %s' % band)
        maps[band] = survey_depth(l, band, etc, nside=nside)

    if savemaps:
        np.savez(output_dir + os.sep + 'depth.npz', 
                 **maps)

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
    

def average_cadence_slow(l, band=None, nside=64, exclude_bands=[], season_gap=100):
    """
    """
    npix = hp.nside2npix(nside)
    cad = np.zeros(npix)
    cad[:] = hp.UNSEEN
    
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


def average_cadence(visits, band=None, nside=64, exclude_bands=[], season_gap=100, compressed=False):
    """
    """
    npix = hp.nside2npix(nside)
    cad_median = np.zeros(npix)
    cad_median[:] = hp.UNSEEN

    cad_mad = np.zeros(npix)
    cad_mad[:] = hp.UNSEEN

    cad_mean = np.zeros(npix)
    cad_mean[:] = hp.UNSEEN

    cad_rms = np.zeros(npix)
    cad_rms[:] = hp.UNSEEN

    cad_gap_10 = np.zeros(npix)
    cad_gap_10[:] = hp.UNSEEN

    cad_gap_15 = np.zeros(npix)
    cad_gap_15[:] = hp.UNSEEN
    
    if band is not None:
        v = visits[visits['band'] == band]
    else:
        v = visits
    
    bins = np.arange(0., season_gap, 1.)
    ipix = np.unique(v['pixel'])
    for j,i in enumerate(ipix):
        if j%100 == 0:
            logging.info('pixel: %d / %d' % (j, len(ipix)))
        idx = v['pixel'] == i
        mjd = v['mjd'][idx]
        dt = mjd[1:] - mjd[:-1]
        #        assert not np.any(dt<=0.)
        dt = dt[dt>0.]
        cad_median[i] = np.median(dt[dt<season_gap])    # np.mean(dt[dt<season_gap])
        cad_mad[i] = mad(dt[dt<season_gap]) # np.std(dt[dt<season_gap])
        cad_mean[i] = np.mean(dt[dt<season_gap])    # np.mean(dt[dt<season_gap])
        cad_rms[i] = np.std(dt[dt<season_gap]) # np.std(dt[dt<season_gap])
        ii = dt>0.
        if ii.sum() <= 1:
            continue
        c, _ = np.histogram(dt[ii], bins=bins)
        if np.any(np.isnan(c)):
            logging.warning(' -> nan detected for pixel %d (ii.sum()=%d, in dt : %d)' % (i, ii.sum(), int(np.any(np.isnan(dt[ii])))))
            continue
        cad_gap_10[i] = np.sum(c[11:]) / float(np.sum(c))
        cad_gap_15[i] = np.sum(c[16:]) / float(np.sum(c))
            
    return cad_mean, cad_rms, cad_median, cad_mad, cad_gap_10, cad_gap_15


def _4_band_map_plot(maps, key, fignum=1, min=0., max=30., nest=1, output_dir=None, plotname=None, 
                     suptitle='', units=''):
    fig = pl.figure(num=1, figsize=(8,6))
    fig.clear()
    pl.suptitle(suptitle)
    for i,band in zip([1,2,3,4],'griz'):
        cad = maps[band + key]
        hp.mollview(cad, min=min, max=max, nest=nest, fig=fignum, sub=(2,2,i), 
                    title='%s [%5.1f %s]' % (band, np.median(cad[cad>0.]), units))
    if output_dir is not None and plotname is not None:
        _savefig(fig, output_dir + os.sep + plotname)


def _1_map_plot(maps, key, fignum=1, min=0., max=30., nest=1, output_dir=None, plotname=None, 
                     suptitle='', units=''):
    fig = pl.figure(num=fignum, figsize=(8,6))
    fig.clear()
    pl.suptitle(suptitle)
    cad = maps[key]
    hp.mollview(cad, min=min, max=max, nest=nest, fig=fignum,
                title='griz [%5.1f %s]' % (np.median(cad[cad>0.]), units))
    if output_dir is not None and plotname is not None:
        _savefig(fig, output_dir + os.sep + plotname)

def plot_average_cadence(data, nside=64, output_dir=None, savemaps=False, exclude_bands=['u', 'y'], filename='cadence.npz', plotname='average_cadence.png'):
    """
    """
    logging.info(' ? visits ? %r' % visits)
    
    # make sure that bands excluded indeed 
    N = len(data)
    idx_excl = np.zeros(N).astype(bool)
    for bn in exclude_bands:
        idx_excl &= (data['band'] == bn)
    if idx_excl.sum() > 0:
        data = data[~idx_excl]    

    # main cadence
    maps = {}
    for i,band in zip([1,2,3,4],'griz'):
        logging.info('processing band: %s' % band)
        cad_mean, cad_rms, cad_median, cad_mad, cad_gap_10, cad_gap_15 = average_cadence(data, band=band, nside=nside) # was visits 
        maps[band] = cad_mean
        maps[band + '_rms'] = cad_rms
        maps[band + '_median'] = cad_median
        maps[band + '_mad'] = cad_mad
        maps[band + '_gap_10'] = cad_gap_10
        maps[band + '_gap_15'] = cad_gap_15

    _4_band_map_plot(maps, key='', fignum=1, output_dir=output_dir, plotname=plotname, 
                     suptitle='mean $\Delta t$ between obs (Full survey)', 
                     units='days')
    _4_band_map_plot(maps, key='_rms', fignum=2, output_dir=output_dir, plotname='rms_' + plotname, 
                     suptitle='std($\Delta t$) between obs (Full survey)', 
                     units='days')
    _4_band_map_plot(maps, key='_gap_10', fignum=3, min=0, max=None, output_dir=output_dir, plotname='gap_10_' + plotname, 
                     suptitle='fraction of 10 day gaps in observervation sequence', 
                     units='')
    _4_band_map_plot(maps, key='_gap_15', fignum=3, min=0, max=None, output_dir=output_dir, plotname='gap_15_' + plotname, 
                     suptitle='fraction of 15 day gaps in observervation sequence', 
                     units='')

        
    # the same, with griz together
    cad_mean, cad_rms, cad_median, cad_mad, cad_gap_10, cad_gap_15 = average_cadence(data, band=None, nside=nside) # was visits
    maps['griz'] = cad_mean
    maps['griz_rms'] = cad_rms
    maps['griz_median'] = cad_median
    maps['griz_mad'] = cad_mad
    maps['griz_gap_10'] = cad_gap_10
    maps['griz_gap_15'] = cad_gap_15

    if savemaps:
        np.savez(output_dir + os.sep + filename,
                 **maps)

    _1_map_plot(maps, key='griz', fignum=1, min=0, max=30, output_dir=output_dir, plotname='all_' + plotname, 
                suptitle='mean $\Delta t$ between obs (Full survey, griz)', units='days')
    _1_map_plot(maps, key='griz_rms', fignum=2, min=0, max=30, output_dir=output_dir, plotname='all_rms_' + plotname, 
                suptitle='std($\Delta t$) between obs (Full survey, griz)', units='days')
    _1_map_plot(maps, key='griz_gap_10', fignum=3, min=0, max=None, output_dir=output_dir, plotname='all_gap_10_' + plotname, 
                suptitle='fraction of $>10$ day gaps in obs sequence (all filters)', units='')
    _1_map_plot(maps, key='griz_gap_15', fignum=4, min=0, max=None, output_dir=output_dir, plotname='all_gap_15_' + plotname, 
                suptitle='fraction of $>15$ day gaps in obs sequence (all filters)', units='')
        


def n_seasons_and_median_duration_slow(l, nside=64, exclude_bands=['u', 'y'], season_gap=100, min_season_length=40.):
    npix = hp.nside2npix(nside)
    n_seasons = np.zeros(npix)
    median_season_duration = np.zeros(npix)
    N = len(l)
    idx_excl = np.zeros(N).astype(bool)
    for bn in exclude_bands:
        idx_excl &= (l['band'] == bn)
    ll = l[~idx_excl]
    
    ipix = np.unique(ll['pixel'])
    mjd = np.floor(ll['mjd'])
    for j,i in enumerate(ipix):
        if j%100 == 0:
            logging.info('pixel: %d / %d' % (j, len(ipix)))
        idx = ll['pixel'] == i
        m = np.unique(mjd[idx])
        dm = m[1:] - m[:-1]
        k = np.where(dm>season_gap)[0]
        season_length = np.array([uu[1:].sum() for uu in np.split(dm, k)])
        season_length = season_length[season_length>min_season_length]
        n_seasons[i] = len(season_length)
        median_season_duration[i] = np.median(season_length)
    return n_seasons, median_season_duration


def n_seasons_and_median_duration(visits, nside=64, exclude_bands=['u', 'y'], season_gap=100, min_season_length=40.):
    """
    .. note: assumes visits are sorted by mjd 
    """
    npix = hp.nside2npix(nside)
    n_seasons = np.zeros(npix)
    median_season_duration = np.zeros(npix)
    
    ipix = np.unique(visits['pixel'])
    for j,i in enumerate(ipix):
        if j%100 == 0:
            logging.info('pixel: %d / %d' % (j, len(ipix)))
        idx = visits['pixel'] == i
        mjd = visits['mjd'][idx]
        dm = mjd[1:] - mjd[:-1]
        k = np.where(dm>season_gap)[0]
        season_length = np.array([uu[1:].sum() for uu in np.split(dm, k)])
        season_length = season_length[season_length>min_season_length]
        n_seasons[i] = len(season_length)
        median_season_duration[i] = np.median(season_length)
    return n_seasons, median_season_duration


def plot_seasons(data, nside=64, exclude_bands=['u', 'y'], season_gap=100, min_season_length=40., output_dir=None, savemaps=True):
    """
    """
    ns, sd = n_seasons_and_median_duration(data, nside=nside,
                                           exclude_bands=exclude_bands,
                                           season_gap=season_gap,
                                           min_season_length=min_season_length)
    ns[ns<=0.] = hp.UNSEEN
    sd[sd<=0.] = hp.UNSEEN
    if savemaps:
        np.savez(output_dir + os.sep + 'seasons.npz', 
                 number_of_seasons=ns, 
                 season_length=sd)
    # number of seasons
    fig = pl.figure(num=1, figsize=(8,6))
    fig.clear()
    #    pl.suptitle('Search seasons (Full survey)')
    # main cadence
    hp.mollview(ns, min=0, nest=1, fig=1, sub=(2,1,1), 
                title='number of seasons [%5.1f]' % (np.median(ns[ns>0.])))
    hp.mollview(sd, min=0, nest=1, fig=1, sub=(2,1,2), 
                title='season length [%5.1f days]' % (np.median(sd[sd>0.])))
    _savefig(fig, output_dir + os.sep + 'seasons.png')


if __name__ == '__main__':
    parser = ArgumentParser(
        description="""global stats for the cadence specified in argumement""")
    parser.add_argument('-O', '--output_dir',
                        default=None,
                        help='output directory')
    parser.add_argument('--nside',
                        default=64, type=int,
                        help='output directory')
    parser.add_argument('obslog', type=str,
                        help='obslog')
    args = parser.parse_args()
    print args
    
    # open the obslog 
    #    cad = np.load(args.cadence_filename)
    with np.load(args.obslog) as f:
        obslog = f['l']
    
    obslog.sort(order='mjd')
    
    logging.info('compressing log info')
    visits = hpixlog.build_visits(obslog, exclude_bands=['u', 'y'])
    visits.sort(order='mjd')    
    logging.info('done')
    
    # global stats
    #    r = cadence_stats(cad, bands='ugrizy', title=args.cadence_filename)

    # numbers of visits
    plot_number_of_visits(obslog, nside=args.nside, max_nv=400, output_dir=args.output_dir, savemaps=True, filename='raw_nvisits.npz', plotname='raw_nvisits.png')
    plot_number_of_visits(visits, nside=args.nside, max_nv=400, output_dir=args.output_dir, savemaps=True, filename='compressed_nvisits.npz', plotname='compressed_nvisits.png')

    # average seeing
    plot_average_seeing(obslog, nside=args.nside, output_dir=args.output_dir, savemaps=True)
    
    # survey depth
    etc = psf.find('LSSTPG')
    plot_survey_depth(obslog, etc, nside=args.nside, maxv=28., output_dir=args.output_dir, savemaps=True)

    # ... aaaaand finally ! (long !)
    if True:
        plot_average_cadence(visits, nside=args.nside, exclude_bands=['u', 'y'], output_dir=args.output_dir, savemaps=True, filename='compressed_cadence.npz', plotname='compressed_cadence.png')
#        plot_average_cadence(obslog, nside=args.nside, exclude_bands=['u', 'y'], output_dir=args.output_dir, savemaps=True, filename='raw_cadence.npz', plotname='raw_cadence.png')
#        plot_seasons(visits, nside=args.nside, output_dir=args.output_dir, savemaps=True)        
        





    # # cadence plot
    # fig = pl.figure(num=1, figsize=(8,6))
    # fig.clear()
    # pl.suptitle('median $\Delta t$ between obs (Full survey)')
    # for i,band in zip([1,2,3,4],'griz'):
    #     cad = maps[band]
    #     hp.mollview(cad, min=0, max=30, nest=1, fig=1, sub=(2,2,i), 
    #                 title='%s [%5.1f days]' % (band, np.median(cad[cad>0.])))
    # _savefig(fig, output_dir + os.sep + plotname)
    
    # # rms plot 
    # fig_rms = pl.figure(num=2, figsize=(8,6))
    # fig_rms.clear()
    # pl.suptitle('MAD($\Delta t$) between obs (Full survey)')
    # for i,band in zip([1,2,3,4],'griz'):
    #     cad_rms = maps[band + '_rms']
    #     hp.mollview(cad_rms, min=0, max=30, nest=1, fig=2, sub=(2,2,i), 
    #                 title='%s [%5.1f days]' % (band, np.median(cad_rms[cad_rms>0.])))
    # _savefig(fig_rms, output_dir + os.sep + 'rms_' + plotname)

    # # gap_10 plot
    # fig_rms = pl.figure(num=2, figsize=(8,6))
    # fig_rms.clear()
    # pl.suptitle('MAD($\Delta t$) between obs (Full survey)')
    # for i,band in zip([1,2,3,4],'griz'):
    #     cad_rms = maps[band + '_rms']
    #     hp.mollview(cad_rms, min=0, max=30, nest=1, fig=2, sub=(2,2,i), 
    #                 title='%s [%5.1f days]' % (band, np.median(cad_rms[cad_rms>0.])))
    # _savefig(fig_rms, output_dir + os.sep + 'rms_' + plotname)
    
    # # hist plot
    # fig_hist = pl.figure(num=3, figsize=(8,6))
    # fig_hist.clear()
    # pl.suptitle("Histogram of $\Delta t$'s")
    # for i,band in zip([1,2,3,4],'griz'):
    #     cad_hist = maps[band + '_hist']
    #     ax = pl.subplot(2,2,i)
    #     pl.plot(cad_hist, 'bo')
    #     pl.title(band)
    #     pl.xlabel('$\delta t$ [days]')
    # _savefig(fig_hist, output_dir + os.sep + 'hist_' + plotname)


    # fig = pl.figure(num=1, figsize=(8,6))
    # fig.clear()
    # pl.suptitle('median $\Delta t$ between obs (Full survey, griz)')

    # hp.mollview(cad, min=0, max=30, nest=1, fig=1,
    #             title='%s [%5.1f days]' % (band, np.median(cad[cad>0.])))
    # _savefig(fig, output_dir + os.sep + 'all_' + plotname)
    
    # fig_rms = pl.figure(num=2, figsize=(8,6))
    # fig_rms.clear()
    # pl.suptitle('MAD($\Delta t$) between obs (Full survey, griz)')
    # hp.mollview(cad_rms, min=0, max=30, nest=1, fig=2,
    #             title='%s [%5.1f days]' % (band, np.median(cad_rms[cad_rms>0.])))
    # _savefig(fig_rms, output_dir + os.sep + 'all_rms_' + plotname)

    # # hist plot
    # fig_hist = pl.figure(num=3, figsize=(8,6))
    # fig_hist.clear()
    # pl.suptitle("Histogram of $\Delta t$'s")
    # pl.plot(hist, 'bo')
    # pl.xlabel('$\delta t$ [days]')
    # _savefig(fig_hist, output_dir + os.sep + 'all_hist_' + plotname)
