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

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.interactive(0)

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

    def _update_cache(self, mjd_min, mjd_max):
        if self.cache is None or \
           (mjd_max > self.cache_mjd_max) or \
           (mjd_min < self.cache_mjd_min):
            d = self.pixel_obs
            mjd_min, mjd_max = mjd_min-20., mjd_max + 50.
            logging.debug('cache fault -> updating [%7.0f,%7.0f]' % \
                          (mjd_min, mjd_max))
            idx = (d['mjd']>=mjd_min) & (d['mjd']<=mjd_max)
            self.cache = d[idx]
            self.cache_mjd_max = mjd_max
        
    def select_window(self, mjd, **kwargs):
        z = kwargs.get('z', 0.)
        mjd_min, mjd_max = mjd + self.rf_phase_range * (1.+z)
        self._update_cache(mjd_min, mjd_max)
        logging.debug('z: %5.2f mjd: [%f,%f]' % (z, mjd_min, mjd_max))
        idx = (self.cache['mjd']>mjd_min) & (self.cache['mjd']<mjd_max)
        if 'band' in kwargs:    
            idx &= (self.cache['band'] == kwargs['band'])
        return self.cache[idx]

    def accept(self, pxobs):
        t = np.array([f(pxobs) for f in self.accept])
        return np.all(t)
    
    def cadence(self, pxobs, **kwargs):
        """
        Note: if there are back-to-back exposures of the same field
        taken during one single night, then they should be counted at
        one single visit.
        """
        z = kwargs.get('z', 0.)
        mjd_min, mjd_max = pxobs['mjd'].min(), pxobs['mjd'].max()
        dt = (mjd_max-mjd_min) / (1.+z)
        c = np.bincount(pxobs['pixel'], minlength=self.npix).astype(float)
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
                    min=kwargs.get('vmin', 0.), 
                    max=kwargs.get('vmax', 1.),
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
                    min=kwargs.get('vmin', 0.), 
                    max=kwargs.get('vmax', 1.),
                    sub=kwargs.get('sub', None),
                    cbar=kwargs.get('cbar', True))
        plt.title(kwargs.get('title', ''))
        if 'dump_plot_dir' in kwargs:
            fig = plt.gcf()
            fig.savefig(kwargs['dump_plot_dir'] + os.sep + '%05d.png' % self.fig_odometer)
            self.fig_odometer += 1
            fig.clear()        
        return plt.gcf()


def movie(l, zmax=0.5, bands="gri", nside=64, dump_plot_dir=None):
    """
    make a movie of the fields that pass the cuts, up to a redshift z.
    """

    m = Metrics(l)
    
    for mjd in np.arange(m.mjd.min(), m.mjd.max()+1):
        zmax = np.zeros(m.npix)
        nsn  = np.zeros(m.npix)
        
        # check that the sampling is ok at z=0
        s = m.select_window(mjd, z=0.)
        c = m.cadence(s, z=0.)
        first, last = m.first_last_visits(mjd, s, z=0.)
        c[c<0.5] = 0.
        c[first==0.] = 0.
        c[last==0.] = 0.
        c0_ok = c>0.
        
        # loop over the redshift range, and check the resolution in
        # color as a function of redshift. Store the highest redshift
        # that passes the cuts 
        for z in np.arange(0.3, 0.51, 0.01)[::-1]:
            s = m.select_window(mjd, z=z)
            cz = m.cadence(s, z=z)
            firstz, lastz = m.first_last_visits(mjd, s, z=z)
            cz[(cz<0.5)] = 0.
            cz[(firstz==0.)] = 0.
            cz[(lastz==0)] = 0.
            cz *= c0_ok
            zmax[(cz>0) & (zmax==0.)] = z
            c[c==0] = cz[c==0]
        # update the number of supernovae for that day 
        # we update (1) a map that contains the total
        # number of SNe and (2) a NTuple that contains
        # mjd, nsn, zmax
        
            
        m.plot_map(first, fig=1, vmin=0., vmax=1.25, sub=221, cbar=False)
        m.plot_map(last, fig=1, vmin=0., vmax=1.25, sub=222, cbar=False)
        m.plot_map(zmax, fig=1, vmin=0., vmax=0.5, sub=223, cbar=True)        
        m.plot_cadence(c, fig=1, dump_plot_dir=dump_plot_dir, title='%6.0f' % mjd,
                       vmin=0.,
                       vmax=1.25,
                       min_cadence=0.5,
                       sub=224,
                       cbar=True)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="make a movie out of a cadence")
    parser.add_argument('-O', '--output-dir',
                        default='./',
                        help='output directory (for the plots)')
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
        
    #    m = Metrics(l)
    movie(l, bands='gri', dump_plot_dir=args.output_dir)
    
    
        

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
            
            