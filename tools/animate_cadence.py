#!/usr/bin/env python

"""
"""

import os
import os.path as op
from exceptions import ValueError
import logging
logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',
                    level=logging.DEBUG)
import argparse
from mx.DateTime import DateTimeFromMJD

import matplotlib
from matplotlib.cm import jet, bwr, PiYG
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.interactive(0)

matplotlib.rcParams['text.latex.preamble'] = [
    r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
    r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
    r'\usepackage{helvet}',    # set the normal font here
    r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
    r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
]

import numpy as np
import healpy as hp

from saunerie import psf, snsim
from pycosmo import cosmo


instrument_name = 'LSSTPG'


def main(filename='h.npy.npz'):
    f = np.load(filename)
    m = f['m']
    hp.mollview(m, nest=1)

    l = f['l']
    return l

maps = [jet, bwr, PiYG]
for m in maps:
    m.set_under('w')
    m.set_over(m(1.0))
    m.set_under('w')
    m.set_bad('gray')
    #bwr.set_under('w')
    #PiYG.set_under('w')



class Metrics(object):
    """
    
    """
    def __init__(self, pxlobs, z=0.5, rf_phase_range=(-20., 45.), nside=64, etc=None, lc_template=None, model_filename=None):
        """
        Constructor. Takes a file 
        """
        self.pixel_obs = pxlobs
        self.mjd = np.unique(np.floor(pxlobs['mjd'])).astype(int)
        self.mjd = self.mjd[self.mjd>0]
        self.nside, self.npix = nside, hp.nside2npix(nside)
        self.rf_phase_range = np.array(rf_phase_range)
        self.min_rf_phase_range = np.array(rf_phase_range) + np.array([10., -10.])
        self.accept = []
        self.fig_odometer = 0
        self.cache_mjd_min = 0.
        self.cache_mjd_max = 0.
        self.cache = None
        
        # add an etc
        logging.info('loading ETC')
        self.etc = psf.find(instrument_name) if etc is None else etc
        suffix = self.etc.instrument.name + '::'
        self.pixel_obs['bandname'][:] = np.core.defchararray.add(suffix, self.pixel_obs['band'])
        
        # add a LC model
        self.sn = np.rec.fromarrays(([0.] * 7),
                                    names=['z', 'dL', 'DayMax', 'X1', 'Color', 'ra', 'dec'])
        self.sn['X1'], self.sn['Color'] = -2., 0.2
        logging.info('instantiating Cosmo model')
        self.cosmo = cosmo.CosmoLambda()
        logging.info('loading SALT2 from: %s' % model_filename)
        #        self.lcmodel = snsim.init_lc_model(np.unique(self.pixel_obs['bandname']),
        #                                           model_filename)
        logging.info('done.')        
        
        # add something to manage the LC templates
        #        print lc_template
        if lc_template is not None:
            self.lc_template_cache = self._build_lc_template_cache(lc_template)
        else:
            self.lc_template_cache = None

    def _build_lc_template_cache(self, lc_template):
        logging.info('building light curve template cache from: %s' % lc_template)
        d = np.load(lc_template) if type(lc_template) is str else lc_template
        ret = {}
        bands = np.unique(d['band'])
        redshifts = np.unique(d['z'])
        for b in bands:
            for z in redshifts:
                key = '%5.3f-%s' % (z,b)
                idx = (d['band'] == b) & (d['z'] == z)
                if idx.sum() == 0:
                    logging.warning('skipped %r (no entries)' % key)
                    continue
                t,v = d[idx]['mjd'], d[idx]['val']
                ret[key] = (t,v)
        logging.info('Done. Cache has %d entries' % len(ret))
        return ret

    def get_lc_template(self, band=None, z=None):
        key = '%5.3f-%s' % (z,band)
        r = self.lc_template_cache.get(key, (None,None))
        #        print key, self.lc_template_cache.keys()
        return r
    
    def _update_cache(self, mjd_min, mjd_max, **kwargs):
        if self.cache is None or \
           (mjd_max > self.cache_mjd_max) or \
           (mjd_min < self.cache_mjd_min):
            d = self.pixel_obs
            mjd_min, mjd_max = mjd_min-20., mjd_max + 50.
            logging.debug('cache fault -> updating [%7.0f,%7.0f]' % \
                          (mjd_min, mjd_max))
            idx = (d['mjd']>=mjd_min) & (d['mjd']<=mjd_max)
            if 'exclude_bands' in kwargs:
                for b in kwargs['exclude_bands']:
                    idx &= (d['band'] != b)
            self.cache = d[idx]
            self.cache['mjd'] = np.floor(self.cache['mjd'])
            self.cache_mjd_max = mjd_max
            visits = np.rec.fromarrays((self.cache['pixel'], self.cache['mjd'], self.cache['band']),
                                       names=['pixel', 'mjd', 'band'])
            self.visit_cache = np.unique(visits)
        
    def select_window(self, mjd, **kwargs):
        z = kwargs.get('z', 0.)
        mjd_min, mjd_max = mjd + self.rf_phase_range * (1.+z)
        self._update_cache(mjd_min, mjd_max, **kwargs)
        logging.debug('z: %5.2f mjd: [%f,%f]' % (z, mjd_min, mjd_max))
        idx = (self.cache['mjd']>mjd_min) & (self.cache['mjd']<mjd_max)
        if 'band' in kwargs:    
            idx &= (self.cache['band'] == kwargs['band'])
        r = self.cache[idx]
        
        idx = (self.visit_cache['mjd']>mjd_min) & (self.visit_cache['mjd']<mjd_max)
        if 'band' in kwargs:    
            idx &= (self.visit_cache['band'] == kwargs['band'])
        u = self.visit_cache[idx]
        
        return r,u

    def accept(self, pxobs):
        t = np.array([f(pxobs) for f in self.accept])
        return np.all(t)
    
    def cadence(self, mjd_px, **kwargs):
        """
        Note: if there are back-to-back exposures of the same field
        taken during one single night, then they should be counted at
        one single visit.
        """
        z = kwargs.get('z', 0.)
        mjd_min, mjd_max = mjd_px['mjd'].min(), mjd_px['mjd'].max()
        dt = (mjd_max-mjd_min) / (1.+z)
        c = np.bincount(mjd_px['pixel'], minlength=self.npix).astype(float)
        #        c = np.bincount(pxobs['pixel'], minlength=self.npix).astype(float)
        c /= dt
        return c

    def first_last_visits(self, mjd, pxobs, **kwargs):
        z = kwargs.get('z', 0.)
        mjd_min, mjd_max = mjd + self.min_rf_phase_range * (1.+z)
        idx_min = pxobs['mjd']<=mjd_min
        idx_max = pxobs['mjd']>=mjd_max
        c_min = np.bincount(pxobs['pixel'][idx_min], minlength=self.npix).astype(float)
        c_max = np.bincount(pxobs['pixel'][idx_max], minlength=self.npix).astype(float)        
        return c_min, c_max

    def amplitude_snr(self, mjd, bandname, z, px_obs):
        c = np.zeros(self.npix)
        p_tpl, v_tpl = self.get_lc_template(band=bandname,z=z)
        if p_tpl is None or v_tpl is None:
            logging.warning('unable to find LC template for (%s,%f)' % (bandname,z))
            return None
        #        self.sn['z'] = z
        #        self.sn['dL'] = self.cosmo.dL(z)
        d = px_obs[px_obs['bandname'] == bandname]
        f5 = self.etc.mag_to_flux(d['m5'], d['bandname'])
        Li = np.interp(d['mjd']-mjd, p_tpl, v_tpl, left=0., right=0.)
        #        Li = self.lcmodel(self.sn, d['mjd'], d['bandname'], jac=False)
        v = (Li/f5)**2
        snr = np.bincount(d['pixel'], weights=v, minlength=self.npix)
        #        print ' --> ', np.median(snr)
        return 5. * np.sqrt(snr)

    def cut_on_amplitude_snr(self, mjd, z, px_obs, snr_cuts):
        r, observed = np.ones(self.npix), np.zeros(self.npix)
        for band in snr_cuts.keys():
            snr = self.amplitude_snr(mjd, band, z, px_obs)
            #            print np.mean(snr), snr_cuts[band], np.sum(snr<snr_cuts[band])
            r[snr<snr_cuts[band]] = 0.
            observed[snr>0.] = 1
        r *= observed
        return r
    
    def sigma_color(self, mjd, px_obs, z):
        """
        Using SALT2 directly
        """
        pass
    
    def plot_map(self, c, **kwargs):
        """
        plot the cadence on the pixellized sphere. 
        """
        cc = c.copy()
        cc[cc==0.] = hp.UNSEEN            
        hp.mollview(cc, nest=1, fig=kwargs.get('fig', None),
                    min=kwargs.get('vmin', None), 
                    max=kwargs.get('vmax', None),
                    sub=kwargs.get('sub', None),
                    cbar=kwargs.get('cbar', True),
                    cmap=jet)
        # hp.graticule()
        plt.title(kwargs.get('title', ''))
        if 'dump_plot_dir' in kwargs:
            fig = plt.gcf()
            prefix = kwargs.get('prefix', '')
            fig.savefig(kwargs['dump_plot_dir'] + os.sep + prefix + '%05d.png' % self.fig_odometer)
            fig.clear()
        return plt.gcf()

    def plot_cadence(self, c, **kwargs):
        """
        plot the cadence on the pixellized sphere. 
        """
        cc = c.copy()
        if 'min_cadence' in kwargs:
            cc[cc<kwargs['min_cadence']] = 0.
        cc[cc==0.] = hp.UNSEEN            
        hp.mollview(cc, nest=1, fig=kwargs.get('fig', None),
                    min=kwargs.get('vmin', None), 
                    max=kwargs.get('vmax', None),
                    sub=kwargs.get('sub', None),
                    cbar=kwargs.get('cbar', True),
                    cmap=jet)
        # hp.graticule()
        plt.title(kwargs.get('title', ''))
        if 'dump_plot_dir' in kwargs:
            fig = plt.gcf()
            fig.savefig(kwargs['dump_plot_dir'] + os.sep + '%05d.png' % self.fig_odometer)
            fig.clear()        
        return plt.gcf()


class CumulativeNumberOfSNe(object):
    def __init__(self, z, nsn):
        self.z = z
        self.nsn = nsn
    def __call__(self, zz):
        return np.interp(zz, self.z, self.nsn, left=0.)
    
def movie(l, zlim=0.5, zstep=0.01, nside=64, dump_plot_dir=None, nsn_func=None,
          bands=['g', 'r', 'i', 'z'],
          exclude_bands=['u', 'y'],
          vmax_nsn=None,
          min_cadence=0.5,
          lc_template=None,
          salt2=None,
          snr_req=True,
          min_cadence_req=True,
          early_late_meas_req=True):
    """
    make a movie of the fields that pass the cuts, up to a redshift z.
    """
    
    m = Metrics(l, model_filename=salt2, lc_template=lc_template, nside=nside)
    nsn_tot = np.zeros(m.npix)
    nsn_inst = np.zeros(m.npix)
    cadence_tot = np.zeros(m.npix)
    cadence_nhits = np.zeros(m.npix)    
    zmax_tot = np.zeros(m.npix)
    zmax_nhits = np.zeros(m.npix)
    tmp_map = np.zeros(m.npix)

    # median values 
    nsn_tot_history = []
    nsn_inst_history = []
    median_cadence_inst_history = []
    zmax_inst_history = []
    
    
    # loop on the survey mjd -- by steps of 1 day
    for mjd in np.arange(m.mjd.min(), m.mjd.max()+1):
        zmax = np.zeros(m.npix)
        nsn  = np.zeros(m.npix)
        
        # check that the sampling is ok at z=0
        s,u = m.select_window(mjd, z=0., bands=bands, exclude_bands=exclude_bands)
        c = m.cadence(u, z=0.)
        if min_cadence_req:
            c[c<min_cadence] = 0.
            
        if early_late_meas_req:
            first, last = m.first_last_visits(mjd, u, z=0.)
            c[first==0.] = 0.
            c[last==0.] = 0.
        c0_ok = c>0.
        
        # loop over the redshift range, and check the resolution in
        # color as a function of redshift. Store the highest redshift
        # that passes the cuts 
        #        for z in np.arange(0.1, zlim+zstep, zstep)[::-1]:
        for z in np.arange(0.05, zlim+zstep, zstep)[::-1]:
            # select the window 
            s,u = m.select_window(mjd, z=z, exclude_bands=exclude_bands)
            
            # average cadence
            # note: explore median dt
            cz = m.cadence(u, z=z)
            
            # cut in cadence 
            if min_cadence_req:
                cz[(cz<min_cadence)] = 0.
            
            # observations before -15 and after +30 ? 
            if early_late_meas_req:
                firstz, lastz = m.first_last_visits(mjd, u, z=z)
                cz[(firstz==0.)] = 0.
                cz[(lastz==0)] = 0.
                
            # cut on the last visit
            cz *= c0_ok

            # cut on sigma amplitude
            snr_g = snr_r = snr_i = snr_z = None
            if np.abs(z-0.3) <= 0.01:
                snr_g = m.amplitude_snr(mjd, instrument_name + '::g', z, s)
                snr_r = m.amplitude_snr(mjd, instrument_name + '::r', z, s)
                snr_i = m.amplitude_snr(mjd, instrument_name + '::i', z, s)
                snr_z = m.amplitude_snr(mjd, instrument_name + '::z', z, s)                
            if z <= 0.3:
                snr_ok = m.cut_on_amplitude_snr(mjd, z, s, 
                                                snr_cuts = {instrument_name + '::g': 30., 
                                                            instrument_name + '::r': 40., 
                                                            instrument_name + '::i': 30., 
                                                            instrument_name + '::z': 20.})
            else:
                snr_ok = m.cut_on_amplitude_snr(mjd, z, s, 
                                                snr_cuts = {instrument_name + '::r': 40., 
                                                            instrument_name + '::i': 30., 
                                                            instrument_name + '::z': 20.})
            
            # update max-z map 
            idx = zmax == 0.
            if min_cadence_req:
                idx &= (cz>0)
            logging.info('SNR ? %r' % snr_req)
            if snr_req:
                logging.info('applying SNR requirement')
                idx &= (snr_ok>0)
                #            zmax[(cz>0) & (snr_ok>0.) & (zmax==0.)] = z
            zmax[idx] = z
            c[c==0] = cz[c==0]
        # update the number of supernovae for that day 
        # we update (1) a map that contains the total
        # number of SNe and (2) a NTuple that contains
        # mjd, nsn, zmax
        if nsn_func is not None:
            nsn_inst[:] = 0.
            nsn_inst[zmax>0.] = nsn_func(zmax[zmax>0])
            nsn_tot[zmax>0.] += nsn_inst[zmax>0.]
        else:
            logging.warning('no function to compute number of SNe')

        # update the cumulative maps
        cadence_tot += c
        cadence_nhits[c>0] += 1
        zmax_tot += zmax
        zmax_nhits[zmax>0] += 1
        
        #        m.plot_map(first, fig=1, vmin=0., vmax=1.25, sub=221, cbar=False)
        #        m.plot_map(last, fig=1, vmin=0., vmax=1.25, sub=222, cbar=False)
        fig = plt.figure(1, figsize=(15.,7.5))
        human_date = DateTimeFromMJD(mjd).strftime('%Y-%m-%d')
        fig.suptitle('[%s  mjd=%6.0f]' % (human_date, mjd))
        m.plot_map(nsn_tot, fig=1, sub=231, vmin=0., vmax=vmax_nsn, cbar=True, title='$N_{SNe}: %6.0f$ (tot)' % nsn_tot.sum())
        nsn_tot_history.append((mjd,nsn_tot.sum()))
        tmp_map[:] = hp.UNSEEN ; idx = zmax_nhits>0
        tmp_map[idx] = zmax_tot[idx] / zmax_nhits[idx]
        med = np.median(tmp_map[tmp_map>0])
        m.plot_map(tmp_map, fig=1, sub=232, vmin=0., vmax=0.5, cbar=True, title='$z_{max}$ (avg) [%4.2f]' % (med if ~np.isnan(med) else 0))
        tmp_map[:] = hp.UNSEEN ; idx = cadence_nhits>0
        tmp_map[idx] = cadence_tot[idx] / cadence_nhits[idx]
        med = np.median(tmp_map[tmp_map>0])
        m.plot_map(tmp_map, fig=1, sub=233, vmin=0., vmax=1., cbar=True, title='cadence [day$^{-1}$] (avg) [%4.2f]' % (med if ~np.isnan(med) else 0))
        
        m.plot_map(nsn_inst, fig=1, sub=234, vmin=0., vmax=0.015, cbar=True, title='$N_{SNe}: %4.0f$' % nsn_inst.sum())
        nsn_inst_history.append((mjd,nsn_inst.sum()))
        med = np.median(zmax[zmax>0])
        m.plot_map(zmax, fig=1, vmin=0., vmax=0.5, sub=235, cbar=True, title='$z_{max}$ [%4.2f]' % (med if ~np.isnan(med) else 0))
        zmax_inst_history.append((mjd,(med if ~np.isnan(med) else 0)))
        med = np.median(c[c>0])
        m.plot_cadence(c, fig=1, dump_plot_dir=dump_plot_dir, 
                       vmin=0.,
                       vmax=1.,
                       min_cadence=min_cadence,
                       sub=236,
                       title='cadence [day$^{-1}$] [%4.2f]' % (med if ~np.isnan(med) else 0.),
                       cbar=True)
        median_cadence_inst_history.append((mjd,(med if ~np.isnan(med) else 0.)))

        # SNR debug plots 
        fig = plt.figure(2)
        fig.suptitle('[%s  mjd=%6.0f]' % (human_date, mjd))
        if snr_g is not None and snr_r is not None and snr_i is not None and snr_z is not None:
            m.plot_map(snr_g, fig=2, sub=221, vmin=0., vmax=30., cbar=True, title='SNR[g]')
            m.plot_map(snr_r, fig=2, sub=222, vmin=0., vmax=40., cbar=True, title='SNR[r]')
            m.plot_map(snr_i, fig=2, sub=223, vmin=0., vmax=30., cbar=True, title='SNR[i]')        
            m.plot_map(snr_z, fig=2, sub=224, vmin=0., vmax=20., cbar=True, title='SNR[z]', dump_plot_dir=dump_plot_dir, prefix='snr')

        # cadence debug plots

        m.fig_odometer += 1


    # dump history
    nsn_tot_history = np.rec.fromrecords(nsn_tot_history, names=['mjd', 'val'])
    nsn_inst_history = np.rec.fromrecords(nsn_inst_history, names=['mjd', 'val'])
    zmax_inst_history = np.rec.fromrecords(zmax_inst_history, names=['mjd', 'val'])
    median_cadence_inst_history = np.rec.fromrecords(median_cadence_inst_history, names=['mjd', 'val'])
    np.save(dump_plot_dir + os.sep + 'nsn_tot_history.npy', nsn_tot_history)
    np.save(dump_plot_dir + os.sep + 'nsn_inst_history.npy', nsn_inst_history)
    np.save(dump_plot_dir + os.sep + 'zmax_inst_history.npy', zmax_inst_history)
    np.save(dump_plot_dir + os.sep + 'median_cadence_inst_history.npy', median_cadence_inst_history)
    np.save(dump_plot_dir + os.sep + 'nsn_tot.npy', nsn_tot)
    np.save(dump_plot_dir + os.sep + 'zmax_tot.npy', zmax_tot)
    np.save(dump_plot_dir + os.sep + 'cadence_tot.npy', cadence_tot)

        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="make a movie out of a cadence")
    parser.add_argument('-O', '--output-dir',
                        default='./',
                        help='output directory (for the plots)')
    parser.add_argument('--nsn',
                        default=None,
                        help='file containing the tabulated cumulative number of SNe')
    parser.add_argument('--zmax', 
                        default=0.5, type=float,
                        help='highest redshift to test')
    parser.add_argument('--nside',
                        default=128, type=int,
                        help='nside to use in the analysis')
    parser.add_argument('--ebv-mask', 
                        default=None, dest='ebv_mask',
                        help='use pixel mask (generally E(B-V) mask)')
    parser.add_argument('--vmax_nsn', 
                        default=None, dest='vmax_nsn', type=float,
                        help='tune the vmax for the nsn_tot map')
    parser.add_argument('--min_cadence', 
                        default=0.5, dest='min_cadence', type=float,
                        help='discard pixels with cadence < min_cadence')
    parser.add_argument('--salt2', 
                        default='salt2.npz', dest='salt2', type=str,
                        help='SALT2 model')
    parser.add_argument('--drop-min-cadence-req', 
                        default=False, action='store_true',
                        help='drop the requirement on cadence')
    parser.add_argument('--drop-snr-req', 
                        default=False, action='store_true',
                        help='drop the requirement on SNR')
    parser.add_argument('--drop-early-late-req', 
                        default=False, action='store_true',
                        help='drop the requirement on having an early and a late measurement')
    parser.add_argument('--lc_template', 
                        dest='lc_template', type=str,
                        help='light curve template (generated from SALT2)')
    parser.add_argument('obs_file',
                        help='observation log')
    args = parser.parse_args()
    print args
    
    f = np.load(args.obs_file)
    #    m = f['m']
    #    hp.mollview(m, nest=1)
    l = f['l']
    
    if not op.isdir(args.output_dir):
        os.makedirs(args.output_dir)

    if args.nsn is not None:
        f = np.load(args.nsn)
        nsn_func = CumulativeNumberOfSNe(f['z'], f['nsn'])
    else:
        nsn_func = None
        
    if args.ebv_mask is not None:
        mask = np.load(args.ebv_mask)
        idx = mask[l['pixel']] == 1
        l = l[idx]
        logging.info('stripping masked pixels: %d -> %d' % (len(idx), len(l)))
        
    #    m = Metrics(l)
    movie(l, zlim=args.zmax, nside=args.nside, 
          dump_plot_dir=args.output_dir, 
          nsn_func=nsn_func, vmax_nsn=args.vmax_nsn, 
          bands=['g', 'r', 'i', 'z'],
          salt2=args.salt2,
          lc_template=args.lc_template,
          min_cadence=args.min_cadence,
          snr_req= not args.drop_snr_req,
          min_cadence_req = not args.drop_min_cadence_req,
          early_late_meas_req = not args.drop_early_late_req)
    
    #    movie(l, bands='gri', dump_plot_dir=args.output_dir)
    
    # def accept(m):
    #     # more than 2 points before max
    #     # more than 7 points after max
    #     # at least one point after +30
    #     # at least one point before -15
    #     # cadence larger than 0.25 day^-1 (restframe)
    #     # amplitude rms > xxx
    #     pass
    
    # m = Metrics(l, accept)
    # for mjd in m.mjd:
    #     mjd_min, mjd_max = (mjd-20.) * (1+z), (mjd+45.) * (1.+z)
    #     for bn in bands:
    #         d = m.select(mjd=mjd, band=bn, z=z)
    #         pix = m.accept(d)
    #         cad = m.cadence(d)
    #         snr = m.mu_snr(d)
    #         col_snr = m.col_snr(d)
            
            
