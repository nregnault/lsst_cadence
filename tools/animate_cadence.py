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




def main(filename='h.npy.npz'):
    f = np.load(filename)
    m = f['m']
    hp.mollview(m, nest=1)

    l = f['l']
    return l


class Metrics(object):
    
    def __init__(self, pxlobs, z=0.5, rf_phase_range=(-20., 35.), nside=64):
        """
        Constructor. Takes a file 
        """
        self.pixel_obs = pxlobs
        self.mjd = np.unique(np.floor(pxlobs['mjd'])).astype(int)
        self.mjd = self.mjd[self.mjd>0]
        self.nside, self.npix = nside, hp.nside2npix(nside)
        self.rf_phase_range = np.array(rf_phase_range)
        self.min_rf_phase_range = np.array(rf_phase_range) + np.array([5., -5.])
        self.accept = []
        self.fig_odometer = 0
        #        self.cache_mjd_window = np.array((-50., +150.))
        self.cache_mjd_min = 0.
        self.cache_mjd_max = 0.
        self.cache = None
        
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
                    cbar=kwargs.get('cbar', True))
        plt.title(kwargs.get('title', ''))
        if 'dump_plot_dir' in kwargs:
            fig = plt.gcf()
            fig.savefig(kwargs['dump_plot_dir'] + os.sep + '%05d.png' % self.fig_odometer)
            self.fig_odometer += 1
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
                    cbar=kwargs.get('cbar', True))
        plt.title(kwargs.get('title', ''))
        if 'dump_plot_dir' in kwargs:
            fig = plt.gcf()
            fig.savefig(kwargs['dump_plot_dir'] + os.sep + '%05d.png' % self.fig_odometer)
            self.fig_odometer += 1
            fig.clear()        
        return plt.gcf()


class CumulativeNumberOfSNe(object):
    def __init__(self, z, nsn):
        self.z = z
        self.nsn = nsn
    def __call__(self, zz):
        return np.interp(zz, self.z, self.nsn, left=0.)
    
def movie(l, zmax=0.5, nside=64, dump_plot_dir=None, nsn_func=None, bands=['g', 'r', 'i', 'z'], exclude_bands=['u', 'y'], vmax_nsn=None, min_cadence=0.5):
    """
    make a movie of the fields that pass the cuts, up to a redshift z.
    """
    
    m = Metrics(l)
    nsn_tot = np.zeros(m.npix)
    nsn_inst = np.zeros(m.npix)
    
    for mjd in np.arange(m.mjd.min(), m.mjd.max()+1):
        zmax = np.zeros(m.npix)
        nsn  = np.zeros(m.npix)
        
        # check that the sampling is ok at z=0
        s,u = m.select_window(mjd, z=0., bands=bands, exclude_bands=exclude_bands)
        c = m.cadence(u, z=0.)
        first, last = m.first_last_visits(mjd, u, z=0.)
        c[c<min_cadence] = 0.
        c[first==0.] = 0.
        c[last==0.] = 0.
        c0_ok = c>0.
        
        # loop over the redshift range, and check the resolution in
        # color as a function of redshift. Store the highest redshift
        # that passes the cuts 
        for z in np.arange(0.3, 0.51, 0.01)[::-1]:
            s,u = m.select_window(mjd, z=z, exclude_bands=exclude_bands)
            cz = m.cadence(u, z=z)
            firstz, lastz = m.first_last_visits(mjd, u, z=z)
            cz[(cz<min_cadence)] = 0.
            cz[(firstz==0.)] = 0.
            cz[(lastz==0)] = 0.
            cz *= c0_ok
            zmax[(cz>0) & (zmax==0.)] = z
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
            
        #        m.plot_map(first, fig=1, vmin=0., vmax=1.25, sub=221, cbar=False)
        #        m.plot_map(last, fig=1, vmin=0., vmax=1.25, sub=222, cbar=False)
        fig = plt.gcf()
        human_date = DateTimeFromMJD(mjd).strftime('%Y-%m-%d')
        fig.suptitle('[%s  mjd=%6.0f]' % (human_date, mjd))
        m.plot_map(nsn_tot, fit=1, sub=221, vmin=0., vmax=vmax_nsn, cbar=True, title='$N_{SNe}: %6.0f$' % nsn_tot.sum())
        m.plot_map(nsn_inst, fit=1, sub=222, vmin=0., cbar=False, title='$N_{SNe}[%s]: %4.0f$' % (human_date, nsn_inst.sum()))
        m.plot_map(zmax, fig=1, vmin=0., vmax=0.5, sub=223, cbar=True, title='$z_{max}$')
        m.plot_cadence(c, fig=1, dump_plot_dir=dump_plot_dir, 
                       vmin=0.,
                       vmax=1.,
                       min_cadence=min_cadence,
                       sub=224,
                       title='cadence [day$^{-1}$]',
                       cbar=True)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="make a movie out of a cadence")
    parser.add_argument('-O', '--output-dir',
                        default='./',
                        help='output directory (for the plots)')
    parser.add_argument('--nsn',
                        default=None,
                        help='file containing the tabulated cumulative number of SNe')
    parser.add_argument('--ebv-mask', 
                        default=None, dest='ebv_mask',
                        help='use pixel mask (generally E(B-V) mask)')
    parser.add_argument('--vmax_nsn', 
                        default=None, dest='vmax_nsn', type=float,
                        help='tune the vmax for the nsn_tot map')
    parser.add_argument('--min_cadence', 
                        default=0.5, dest='min_cadence', type=float,
                        help='discard pixels with cadence < min_cadence')
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
    movie(l, dump_plot_dir=args.output_dir, nsn_func=nsn_func, vmax_nsn=args.vmax_nsn, bands=['g', 'r', 'i', 'z'], min_cadence=args.min_cadence)
            

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
            
            
