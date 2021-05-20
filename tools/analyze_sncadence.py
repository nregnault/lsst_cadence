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




        
class TabulatedLcModel(object):
    
    def __init__(self, data, sntype='faint', max_restframe_phase_range=(-20,45)):
        self.data = data
        self.sntype = sntype
        self.cache = self._build_cache(data, sntype)
        self.min_phase = max_restframe_phase_range[0]
        self.max_phase = max_restframe_phase_range[1]

    def _build_cache(self, data, sntype):
        """
        """
        logging.info('building light curve template cache')
        d = np.load(data) if type(data) is str else data
        ret = {}
        bands = np.unique(d['band'])
        redshifts = np.unique(d['z'])
        for b in bands:
            for z in redshifts:
                key = '%5.3f-%s-%s' % (z,b,sntype)
                idx = (d['band'] == b) & (d['z'] == z) & (d['sntype'] == sntype)
                if idx.sum() == 0:
                    logging.warning('skipped %r (no entries)' % key)
                    continue
                t,v = d[idx]['mjd'], d[idx]['val']
                ret[key] = (t,v)
        logging.info('Done. Cache has %d entries' % len(ret))
        return ret
    
    def __call__(self, mjd, peak_mjd, z, band):
        """
        """
        phase = (mjd-peak_mjd) / (1.+z)
        key = '%5.3f-%s-%s' % (z, band, self.sntype)
        x,y = self.cache.get(key, (None, None))
        r = np.interp(phase, x, y, left=0., right=0.)
        r[phase<self.min_phase] = 0.
        r[phase>self.max_phase] = 0.
        return r
        




class Metrics(object):
    
    def __init__(self, sn_rate, nside=64, ):
        self.sn_rate = sn_rate
        self.nside = nside
        self.nsn_inst = None
        self.nsn_tot = None
        self.cadence_inst = None
        self.cadence_tot = None
        self.maxz_inst = None
        self.maxz_tot = None

        self.dNdz = None

    def __call__(self, lc):
        pass



class DetectedTransient(Metrics):
    """Compute number of transients with at least 2 5sigma detections
    """
    
    def __init__(self, **kwargs):
        super(Metrics, self).__init__(**kwargs)        

    def __call__(self, lc):
        pass

class TwoColorsOnTheRise(Metrics):
    """Number of transients with at rise 2 colors on the rise (measured at 0.1 mag)
    """
    def __init__(self, **kwargs):
        super(Metrics, self).__init__(**kwargs)
        pass

    def __call__(self, lc):
        pass


class HDGradeSN(Metrics):
    """Hubble diagram-grade transients
    
    We demand: 
     - at least one visit before -10
     - at least one visit after  +35
     - a minimum cadence of 0.25
     - snr[g]>, snr[r,i,z]> 
    """
    
    def __init__(self, **kwargs):
        super(Metrics, self).__init__(**kwargs)
        self.min_cadence = kwargs.get('min_cadence', 0.25)
        self.snr_cur = kwargs.get('snr_cut', {'g': 30., 'r': 40., 'i': 30., 'z': 20.})
        N = hp.nside2npix(self.nside)
        
        # number of SNe 
        self.nsn_inst = np.zeros(N)
        self.nsn_tot = np.zeros(N)
        
        # cadence 
        self.cadence_inst = np.zeros(N)
        self.cadence_tot = np.zeros(N)
        self.cadence_counts = np.zeros(N)
        
        # zmax
        self.zmax_inst = np.zeros(N)
        self.zmax_tot = np.zeros(N)
        self.zmax_counts = np.zeros(N)
        
    def cut_on_cadence(self, z, lcw):
        """
        """
        if hasattr(self, 'c0'):
            return self.c0
        c = lcw.cadence(z=0.)
        first, last = lcw.first_last_visits(z=z)
        c[c<self.min_cadence] = 0.
        c[first==0.] = 0.
        c[last==0.] = 0.
        if z == 0.:
            self.c0 = c

    def cut_on_snr(self, z, lcw, lcmodel):
        tested_bands = "griz" if z<0.3 else "riz"
        snr_ok = np.ones(self.N)
        for b in tested_bands:
            s = lcw.amplitude_snr(z, b, lcmodel)
            snr_ok[s<self.snr_cut[b]] = 0.
        return snr_ok
            
    def __call__(self, z, lcw, lcmodel):
        """
        """
        # cut on the first / last visits and on the cadence quality 
        c0 = self.cut_on_cadence(0., lcw)
        cz = self.cut_on_cadence(z, lcw)
        
        # cut on the integrated LC SNR
        snr_ok = self.cut_on_snr(z, lcw, lcmodel)

        # update the cadence maps
        
        
        
        # update the zmax map (faint)
        zmax = self.zmax_inst
        zmax[(cz>0.) & (snr_ok>0.) & (zmax==0.)] = z
        self.zmax_tot += zmax
        self.zmax_counts[zmax>0] += 1.

        # update the zmax map (

        
        
        

class HDGradeSNWithTwoColorsOnTheRise(Metrics):
    """HD-grade transients with at least 2 colors on the rise
    """
    def __init__(self, **kwargs):
        super(Metrics, self).__init__(**kwargs)
        pass

    def __call__(self, lc):
        pass
    


# maybe these should be just functions 
class Plotter(object):
    
    def __init__(self):
        pass

    def plot_movie_nsn_map(self, metrics):
        pass

    def plot_final_map(self, metrics):
        pass

    



class CumulativeNumberOfSNe(object):
    def __init__(self, z, nsn):
        self.z = z
        self.nsn = nsn
    def __call__(self, zz):
        return np.interp(zz, self.z, self.nsn, left=0.)

    

        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="make a movie out of a cadence")
    parser.add_argument('-O', '--output-dir',
                        default='./',
                        help='output directory (for the plots)')
    parser.add_argument('--nsn',
                        default=None,
                        help='file containing the tabulated cumulative number of SNe')
    parser.add_argument('--nside',
                        default=128, type=int,
                        help='nside to use in the analysis')
    parser.add_argument('--ebv-mask', 
                        default=None, dest='ebv_mask',
                        help='use pixel mask (generally E(B-V) mask)')
    parser.add_argument('-s', '--timestep', 
                        default=1., type=float, 
                        help='simulation time step in (observer frame) days')
    parser.add_argument('--vmax_nsn', 
                        default=None, dest='vmax_nsn', type=float,
                        help='tune the vmax for the nsn_tot map')
    parser.add_argument('--min_cadence', 
                        default=0.5, dest='min_cadence', type=float,
                        help='discard pixels with cadence < min_cadence')
    parser.add_argument('--salt2', 
                        default='salt2.npz', dest='salt2', type=str,
                        help='SALT2 model')
    parser.add_argument('--lc_template', 
                        dest='lc_template', type=str,
                        help='light curve template (generated from SALT2)')
    parser.add_argument('obs_file',
                        help='observation log')
    args = parser.parse_args()
    print args
    
    # load the pixel observation log 
    logging.info('loading the pixel observing log from %s' % args.obs_file)
    f = np.load(args.obs_file)
    l = f['l']
    
    # get rid of the masked pixels 
    if args.ebv_mask is not None:
        mask = np.load(args.ebv_mask)
        idx = mask[l['pixel']] == 1
        l = l[idx]
        logging.info('stripping masked pixels: %d -> %d' % (len(idx), len(l)))
        
    # compress the pixel log into a single night visit log 
    # - this should be done in the previous segment of the pipeline
    #    visits = compress(l)
    
    # separate the DDF's from the rest of the survey 
    wfd_data, ddf_data, _ = extract_ddf_data(l, ddf_id=1)
    #    del l
    
    # load the light curve template
    logging.info('loading LC templates from: %s' % args.lc_template)
    tpl = np.load(args.lc_template)
    lcmodel_faint = TabulatedLcModel(tpl, sntype='faint')
    lcmodel_normal = TabulatedLcModel(tpl, sntype='normal')

    # output directory 
    if not op.isdir(args.output_dir):
        os.makedirs(args.output_dir)
        
    ####################################################################
    #                                                                  #
    # main loop on the data                                            #
    #                                                                  #
    ####################################################################
    
    # DDF analysis 
    #    ddf_pxlog = PixelObslog(ddf_data, timestep=args.timestep)
    #    metrics = [DetectedTransient(), TwoColorsOnTheRise(), HDGradeSN(), HDGradeSNWithTwoColorsOnTheRise()]
    #    for lcw in ddf_pxlog:
    #        for z in np.arange(0.05, 1.2, 0.1)[::-1]:
    #            for m in metrics:
    #                r = m(lcw)
    #        for m in metrics:
    #            m.dump(args.output_dir)

    # WFD analysis 
    wfd_pxlog = PixelObslog(wfd_data, timestep=args.timestep)
    metrics = [DetectedTransient(), TwoColorsOnTheRise(), HDGradeSN(), HDGradeSNWithTwoColorsOnTheRise()]
    for lcw in wfd_pxlog:
        for z in np.arange(0.05, 0.6, 0.05)[::-1]:
            lcw.compute(z, lcmodel_faint)
            lcw.compute(z, lcmodel_normal)
            #            print lcw.amplitude_snr.keys()
            for m in metrics:
                r = m(lcw)
                    
        # dump the movie frames -- if requested
#        if args.nsn_movie:
#            dump_nsn_movie_frame(m, args.nsn_movie, args.cadence_name)
#        if args.snr_movie:
#            dump_snr_movie_frame(m, args.snr_movie, args.cadence_name)
#        for m in metrics:
#            m.dump(args.output_dir)

    # don't know yet if this movie has any interest
#    if args.dndz_movie:
#        for m_wfd, m_ddf in zip(wfd_metrics, ddf_metrics):
#            dump_dndz_movie_frame(m_wfd, m_ddf, args.dndz_movie,
#                                  args.cadence_name)    


    
            
            





# class Metrics(object):
#     """
    
#     """
#     def __init__(self, pxlobs, z=0.5, rf_phase_range=(-20., 45.), nside=64, etc=None, lc_template=None, model_filename=None):
#         """
#         Constructor. Takes a file 
#         """
#         self.pixel_obs = pxlobs
#         self.mjd = np.unique(np.floor(pxlobs['mjd'])).astype(int)
#         self.mjd = self.mjd[self.mjd>0]
#         self.nside, self.npix = nside, hp.nside2npix(nside)
#         self.rf_phase_range = np.array(rf_phase_range)
#         self.min_rf_phase_range = np.array(rf_phase_range) + np.array([10., -10.])
#         self.accept = []
#         self.fig_odometer = 0
#         self.cache_min_mjd = 0.
#         self.cache_max_mjd = 0.
#         self.cache = None
        
#         # add an etc
#         logging.info('loading ETC')
#         self.etc = psf.find(instrument_name) if etc is None else etc
#         suffix = self.etc.instrument.name + '::'
#         self.pixel_obs['bandname'][:] = np.core.defchararray.add(suffix, self.pixel_obs['band'])
        
#         # add a LC model
#         self.sn = np.rec.fromarrays(([0.] * 7),
#                                     names=['z', 'dL', 'DayMax', 'X1', 'Color', 'ra', 'dec'])
#         self.sn['X1'], self.sn['Color'] = -2., 0.2
#         logging.info('instantiating Cosmo model')
#         self.cosmo = cosmo.CosmoLambda()
#         logging.info('loading SALT2 from: %s' % model_filename)
#         #        self.lcmodel = snsim.init_lc_model(np.unique(self.pixel_obs['bandname']),
#         #                                           model_filename)
#         logging.info('done.')        
        
#         # add something to manage the LC templates
#         #        print lc_template
#         if lc_template is not None:
#             self.lc_template_cache = self._build_lc_template_cache(lc_template)
#         else:
#             self.lc_template_cache = None

#     def _build_lc_template_cache(self, lc_template):
#         logging.info('building light curve template cache from: %s' % lc_template)
#         d = np.load(lc_template) if type(lc_template) is str else lc_template
#         ret = {}
#         bands = np.unique(d['band'])
#         redshifts = np.unique(d['z'])
#         for b in bands:
#             for z in redshifts:
#                 key = '%5.3f-%s' % (z,b)
#                 idx = (d['band'] == b) & (d['z'] == z)
#                 if idx.sum() == 0:
#                     logging.warning('skipped %r (no entries)' % key)
#                     continue
#                 t,v = d[idx]['mjd'], d[idx]['val']
#                 ret[key] = (t,v)
#         logging.info('Done. Cache has %d entries' % len(ret))
#         return ret

#     def get_lc_template(self, band=None, z=None):
#         key = '%5.3f-%s' % (z,band)
#         r = self.lc_template_cache.get(key, (None,None))
#         #        print key, self.lc_template_cache.keys()
#         return r
    
#     def _update_cache(self, mjd_min, mjd_max, **kwargs):
#         if self.cache is None or \
#            (mjd_max > self.cache_max_mjd) or \
#            (mjd_min < self.cache_min_mjd):
#             d = self.pixel_obs
#             mjd_min, mjd_max = mjd_min-20., mjd_max + 50.
#             logging.debug('cache fault -> updating [%7.0f,%7.0f]' % \
#                           (mjd_min, mjd_max))
#             idx = (d['mjd']>=mjd_min) & (d['mjd']<=mjd_max)
#             if 'exclude_bands' in kwargs:
#                 for b in kwargs['exclude_bands']:
#                     idx &= (d['band'] != b)
#             self.cache = d[idx]
#             self.cache['mjd'] = np.floor(self.cache['mjd'])
#             visits = np.rec.fromarrays((self.cache['pixel'], self.cache['mjd'], self.cache['band']),
#                                        names=['pixel', 'mjd', 'band'])
#             self.visit_cache = np.unique(visits)
        
#     def select_window(self, mjd, **kwargs):
#         z = kwargs.get('z', 0.)
#         mjd_min, mjd_max = mjd + self.rf_phase_range * (1.+z)
#         self._update_cache(mjd_min, mjd_max, **kwargs)
#         logging.debug('z: %5.2f mjd: [%f,%f]' % (z, mjd_min, mjd_max))
#         idx = (self.cache['mjd']>mjd_min) & (self.cache['mjd']<mjd_max)
#         if 'band' in kwargs:    
#             idx &= (self.cache['band'] == kwargs['band'])
#         r = self.cache[idx]
        
#         idx = (self.visit_cache['mjd']>mjd_min) & (self.visit_cache['mjd']<mjd_max)
#         if 'band' in kwargs:    
#             idx &= (self.visit_cache['band'] == kwargs['band'])
#         u = self.visit_cache[idx]
        
#         return r,u

#     def accept(self, pxobs):
#         t = np.array([f(pxobs) for f in self.accept])
#         return np.all(t)
    
#     def cadence(self, mjd_px, **kwargs):
#         """
#         Note: if there are back-to-back exposures of the same field
#         taken during one single night, then they should be counted at
#         one single visit.
#         """
#         z = kwargs.get('z', 0.)
#         mjd_min, mjd_max = mjd_px['mjd'].min(), mjd_px['mjd'].max()
#         dt = (mjd_max-mjd_min) / (1.+z)
#         c = np.bincount(mjd_px['pixel'], minlength=self.npix).astype(float)
#         #        c = np.bincount(pxobs['pixel'], minlength=self.npix).astype(float)
#         c /= dt
#         return c

#     def first_last_visits(self, mjd, pxobs, **kwargs):
#         z = kwargs.get('z', 0.)
#         mjd_min, mjd_max = mjd + self.min_rf_phase_range * (1.+z)
#         idx_min = pxobs['mjd']<=mjd_min
#         idx_max = pxobs['mjd']>=mjd_max
#         c_min = np.bincount(pxobs['pixel'][idx_min], minlength=self.npix).astype(float)
#         c_max = np.bincount(pxobs['pixel'][idx_max], minlength=self.npix).astype(float)        
#         return c_min, c_max

#     def amplitude_snr(self, mjd, bandname, z, px_obs):
#         c = np.zeros(self.npix)
#         p_tpl, v_tpl = self.get_lc_template(band=bandname,z=z)
#         if p_tpl is None or v_tpl is None:
#             logging.warning('unable to find LC template for (%s,%f)' % (band,z))
#             return None
#         #        self.sn['z'] = z
#         #        self.sn['dL'] = self.cosmo.dL(z)
#         d = px_obs[px_obs['bandname'] == bandname]
#         f5 = self.etc.mag_to_flux(d['m5'], d['bandname'])
#         Li = np.interp(d['mjd']-mjd, p_tpl, v_tpl, left=0., right=0.)
#         #        Li = self.lcmodel(self.sn, d['mjd'], d['bandname'], jac=False)
#         v = (Li/f5)**2
#         snr = np.bincount(d['pixel'], weights=v, minlength=self.npix)
#         #        print ' --> ', np.median(snr)
#         return 5. * np.sqrt(snr)

#     def cut_on_amplitude_snr(self, mjd, z, px_obs, snr_cuts):
#         r, observed = np.ones(self.npix), np.zeros(self.npix)
#         for band in snr_cuts.keys():
#             snr = self.amplitude_snr(mjd, band, z, px_obs)
#             #            print np.mean(snr), snr_cuts[band], np.sum(snr<snr_cuts[band])
#             r[snr<snr_cuts[band]] = 0.
#             observed[snr>0.] = 1
#         r *= observed
#         return r
    
#     def sigma_color(self, mjd, px_obs, z):
#         """
#         Using SALT2 directly
#         """
#         pass
    
#     def plot_map(self, c, **kwargs):
#         """
#         plot the cadence on the pixellized sphere. 
#         """
#         cc = c.copy()
#         cc[cc==0.] = hp.UNSEEN            
#         hp.mollview(cc, nest=1, fig=kwargs.get('fig', None),
#                     min=kwargs.get('vmin', None), 
#                     max=kwargs.get('vmax', None),
#                     sub=kwargs.get('sub', None),
#                     cbar=kwargs.get('cbar', True),
#                     cmap=jet)
#         # hp.graticule()
#         plt.title(kwargs.get('title', ''))
#         if 'dump_plot_dir' in kwargs:
#             fig = plt.gcf()
#             prefix = kwargs.get('prefix', '')
#             fig.savefig(kwargs['dump_plot_dir'] + os.sep + prefix + '%05d.png' % self.fig_odometer)
#             fig.clear()
#         return plt.gcf()

#     def plot_cadence(self, c, **kwargs):
#         """
#         plot the cadence on the pixellized sphere. 
#         """
#         cc = c.copy()
#         if 'min_cadence' in kwargs:
#             cc[cc<kwargs['min_cadence']] = 0.
#         cc[cc==0.] = hp.UNSEEN            
#         hp.mollview(cc, nest=1, fig=kwargs.get('fig', None),
#                     min=kwargs.get('vmin', None), 
#                     max=kwargs.get('vmax', None),
#                     sub=kwargs.get('sub', None),
#                     cbar=kwargs.get('cbar', True),
#                     cmap=jet)
#         # hp.graticule()
#         plt.title(kwargs.get('title', ''))
#         if 'dump_plot_dir' in kwargs:
#             fig = plt.gcf()
#             fig.savefig(kwargs['dump_plot_dir'] + os.sep + '%05d.png' % self.fig_odometer)
#             fig.clear()        
#         return plt.gcf()







def movie(l, zlim=0.7, zstep=0.01, nside=64, dump_plot_dir=None, nsn_func=None,
          bands=['g', 'r', 'i', 'z'],
          exclude_bands=['u', 'y'],
          vmax_nsn=None,
          min_cadence=0.5,
          lc_template=None,
          salt2=None):
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
    
    # p = Plotter()
    # for block in pxlog:
    #     for m,acc in metrics:
    #         r = m(block)
    #         a.accumulate(r)
    #         if a.do_plot:
    #             p.plot_maps(a)
                
    
    
    # loop on the survey mjd -- by steps of 1 day
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
        for z in np.arange(0.1, zlim+zstep, zstep)[::-1]:
            # select the window 
            s,u = m.select_window(mjd, z=z, exclude_bands=exclude_bands)
            
            # average cadence
            # note: explore median dt
            cz = m.cadence(u, z=z)
            
            # observations before -15 and after +30 ? 
            firstz, lastz = m.first_last_visits(mjd, u, z=z)
            # cut in cadence 
            cz[(cz<min_cadence)] = 0.
            
            # cut on the last visit
            cz[(firstz==0.)] = 0.
            cz[(lastz==0)] = 0.
            cz *= c0_ok

            # cut on sigma amplitude
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
            zmax[(cz>0) & (snr_ok>0.) & (zmax==0.)] = z
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




#    if args.nsn is not None:
#        f = np.load(args.nsn)
#        nsn_func = CumulativeNumberOfSNe(f['z'], f['nsn'])
#    else:
#        nsn_func = None
        
        
#    #    m = Metrics(l)
#    movie(l, nside=args.nside, dump_plot_dir=args.output_dir, nsn_func=nsn_func, vmax_nsn=args.vmax_nsn, bands=['g', 'r', 'i', 'z'],
#          salt2=args.salt2,
#          lc_template=args.lc_template,
#          min_cadence=args.min_cadence)




# def extract_ddf_data(pxlobs, ddf_id=5, take_all_pixels=False):
#     """Separate the DDF fields (healpixs) from the WFD fields (healpixs)
#     """
#     logging.info('separating the DDF fields (healpix) from the WFD/NGS/SGC')
#     idx = pxlobs['proposal'] == ddf_id
#     ddf_pixels = np.unique(pxlobs['pixel'][idx])    
#     if not take_all_pixels:
#         ddf_data = pxlobs[idx]
#         wfd_data = pxlobs[~idx]
#     else:
#         logging.info('DDF pixels...')
#         logging.info('%d DDF pixels identified' % len(ddf_pixels))
#         idx = np.in1d(pxlobs['pixel'], ddf_pixels)
#         ddf_data = pxlobs[idx]
#         wfd_data = pxlobs[~idx]
#     logging.info('WFD/NGS/SGC: %d obs, DDF: %d obs' % (len(wfd_data), len(ddf_data)))
#     return wfd_data, ddf_data, ddf_pixels


# def build_visits(data, exclude_bands=None):
#     """Ligher index, containing pixel/band/night(mjd) information.
    
#     Build a lighter structure, that contains just the information
#     needed to count the number of individual nights during which a
#     pixel was observed.
#     """
#     logging.info('build night/band array...')
#     mjd = np.floor(data['mjd']).astype(int)
#     visits = np.rec.fromarrays((data['pixel'].astype(int), mjd, data['band']),
#                                names=['pixel', 'mjd', 'band'])
#     logging.info('unique...')
#     visits = np.unique(visits)
#     logging.info('done...')
    
#     idx = np.zeros(len(visits)).astype(bool)
#     if exclude_bands is not None:
#         for b in exclude_bands:
#             idx |= (visits['band'] == b)
#     logging.info('%d visits selected in bands: %r' % ((~idx).sum(), np.unique(visits[~idx]['band'])))
#     return visits[~idx]


# class PixelLightcurveWindow(object):
#     """Manages the fraction around 
    
#     """
    
#     def __init__(self, peak_mjd, pixel_obs, f5, visits,
#                  map_npix=hp.nside2npix(64),
#                  min_restframe_phase_range=(-10, +30.),
#                  max_restframe_phase_range=(-20, +45.), 
#                  etc=None, unique_bandnames=None):
#         self.peak_mjd = peak_mjd
#         self.f5 = f5
#         self.pixel_obs = pixel_obs
#         self.visits = visits
#         self.npix = map_npix
#         self.min_restframe_phase_range = np.array(min_restframe_phase_range)
#         self.max_restframe_phase_range = np.array(max_restframe_phase_range)
#         self.etc = etc
#         self.unique_bandnames = unique_bandnames

#     def __len__(self):
#         return len(self.visits)

    
#     def compute(self, z, lcmodel):
#         self.z = z
#         self.snr, self.amplitude_snr = {}, {}
#         for b in self.unique_bandnames:
#             key = b + '-' + lcmodel.sntype
#             self.amplitude_snr[key] = self.compute_amplitude_snr(z, b, lcmodel)
#         key = lcmodel.sntype
#         self.snr[key] = self.compute_snr(z, lcmodel)

#     def first_last_visits(self, z=0.):
#         """visits before and after a minimum restframe phase range. 
#         """
#         mjd_min_up,  mjd_max_low = self.peak_mjd + self.min_restframe_phase_range * (1.+z)
#         mjd_min_low, mjd_max_up  = self.peak_mjd + self.max_restframe_phase_range * (1.+z)
#         idx_min = (self.visits['mjd']>=mjd_min_low) & (self.visits['mjd']<=mjd_min_up)
#         idx_max = (self.visits['mjd']>=mjd_max_low) & (self.visits['mjd']<=mjd_max_up)
        
#         c_min = np.bincount(self.visits['pixel'][idx_min], minlength=self.npix).astype(float)
#         c_max = np.bincount(self.visits['pixel'][idx_max], minlength=self.npix).astype(float)        
#         return c_min, c_max

#     def restframe_cadence(self, z, band=None):
#         """observation cadence in the SN restframe
#         """
#         mjd_min, mjd_max = self.visits['mjd'].min(), self.visits['mjd'].max()
#         dt = (mjd_max-mjd_min) / (1.+z)
#         if band is None:
#             c = np.bincount(self.visits['pixel'], minlength=self.npix).astype(float)
#         else:
#             idx = self.visits['band'] == band
#             c = np.bincount(self.visits['pixel'][idx], minlength=self.npix).astype(float)
#         c /= dt
#         return c
        
#     def compute_snr(self, z, lcmodel):
#         """given a LC model and a redshift return the SNR of each point 
#         """
#         r = np.zeros(len(self.pixel_obs))
#         d, f5 = self.pixel_obs, self.f5
#         for b in self.unique_bandnames:
#             idx = self.pixel_obs['bandname'] == b
#             Li = lcmodel(d['mjd'][idx], self.peak_mjd, z, b)
#             r[idx] = 5. * Li / f5[idx]
        
#         return r
            
#     def compute_amplitude_snr(self, z, bandname, lcmodel):
#         """
#         """
#         c = np.zeros(self.npix)
#         idx = self.pixel_obs['bandname'] == bandname
#         d, f5 = self.pixel_obs[idx], self.f5[idx]
#         Li = lcmodel(d['mjd'], self.peak_mjd, z, bandname)
#         v = (Li/f5)**2
#         snr = np.bincount(d['pixel'], weights=v, minlength=self.npix)
#         return 5. * np.sqrt(snr)
        
#     def compute_sigma_color(self, z, lcmodel):
#         pass


# class PixelObslog(object):
#     """Iterable interface with a Pixel observation log 


#     .. note : Since the typical pixel observation logs we deal with
#           are large (10 years of survey => 2,4M exposures, nside=64 =>
#           24,000 pixels in the southern sky, resulting in about 25M
#           entries in the pixel log), this class manages a cache, in
#           order to accelerate the requests.
#     """
    
#     def __init__(self, pxlobs, timestep=1, restframe_phase_range=(-25., 45.), 
#                  max_redshift=1.2, exclude_bands=['u', 'y'], instrument_name='LSSTPG'):
#         """
#         """
#         self.pixel_obs = self._clean_log(pxlobs,
#                                          instrument_name=instrument_name,
#                                          exclude_bands=exclude_bands)
        
#         self.visits = build_visits(self.pixel_obs, exclude_bands=exclude_bands)
#         self.mjd = np.unique(np.floor(pxlobs['mjd'])).astype(int)
#         self.mjd = self.mjd[self.mjd>0] # why that ? 
#         self.min_mjd, self.max_mjd = np.min(self.mjd), np.max(self.mjd)
#         self.timestep = timestep
#         self.cache_margin = np.array((-30., 100.))
        
#         logging.info('loading ETC')
#         self.etc = psf.find(instrument_name)
#         self.f5 = self.etc.mag_to_flux(self.pixel_obs['m5'], self.pixel_obs['bandname'])
        
#         self.restframe_phase_range = np.array(restframe_phase_range)
#         self.sliding_window = np.array(restframe_phase_range) * (1. + max_redshift)
#         self.max_redshift = max_redshift
#         self.exclude_bands = exclude_bands
#         logging.info('band names... ')
#         self.unique_bandnames = np.unique(self.pixel_obs['bandname'])
#         self.reinit()
#         logging.info('ready to go')

#     def _clean_log(self, data, instrument_name, exclude_bands):
#         """
#         """
#         idx = np.zeros(len(data)).astype(bool)
#         for b in exclude_bands:
#             idx |= (data['band'] == b)
#         logging.info('select: %d measurements out of %d' % ((~idx).sum(), len(idx)))
#         d = data[~idx]
        
#         suffix = instrument_name + '::'
#         d['bandname'][:] = np.core.defchararray.add(suffix, d['band'])
#         return d
        
#     def reinit(self):
#         """re-initialize the iterator to zero.
#         """
#         self.fig_odometer = 0
#         self.cache_min_mjd = 0.
#         self.cache_max_mjd = 0.
#         self.cache = None
#         self.current_mjd = self.min_mjd

#     def __iter__(self):
#         """iterator function. 
#         """
#         self.reinit()
#         return self 

#     def next(self):
#         """the next function in the iterator 

#         .. note : to be renamed __next__ when upgrading to python 3
#         """
#         dt = self.current_mjd + self.sliding_window
#         self._update_cache(dt[0], dt[1], exclude_bands=self.exclude_bands)
#         self.current_mjd += self.timestep
#         r = self.select_window(self.current_mjd)
#         if len(r) == 0 and (self.current_mjd > self.max_mjd):
#             raise StopIteration
#         return r
        
#     # def _build_visits(self, exclude_bands=None):
#     #     """Ligher index, containing pixel/band/night(mjd) information.
        
#     #     Build a lighter structure, that contains just the information
#     #     needed to count the number of individual nights during which a
#     #     pixel was observed.
#     #     """
#     #     logging.info('build night/band array...')
#     #     mjd = np.floor(self.pixel_obs['mjd']).astype(int)
#     #     visits = np.rec.fromarrays((self.pixel_obs['pixel'].astype(int), mjd, self.pixel_obs['band']),
#     #                                names=['pixel', 'mjd', 'band'])
#     #     logging.info('unique...')
#     #     visits = np.unique(visits)
#     #     logging.info('done...')
        
#     #     idx = np.zeros(len(visits)).astype(bool)
#     #     if exclude_bands is not None:
#     #         for b in exclude_bands:
#     #             idx |= (visits['band'] == b)
#     #     logging.info('%d visits selected in bands: %r' % ((~idx).sum(), np.unique(visits[~idx]['band'])))
#     #     return visits[~idx]
                
#     def _update_cache(self, min_mjd, max_mjd, **kwargs):
#         """Update the internal pixel log cache

#         Since the pixel log is pretty large, we keep a smaller
#         fraction of it, that is enough for ~50 iterations. This speeds
#         up significantly the searches.
        
#         """
#         if self.cache is not None and \
#            (max_mjd <= self.cache_max_mjd) and \
#            (min_mjd >= self.cache_min_mjd):
#             return 
#         d = self.pixel_obs
#         min_mjd = min_mjd + self.cache_margin[0]
#         max_mjd = max_mjd + self.cache_margin[1]
#         logging.debug('cache fault -> updating [%7.0f,%7.0f]' % \
#                       (min_mjd, max_mjd))
        
#         # select observations + cut on unused bands
#         idx = (d['mjd']>=min_mjd) & (d['mjd']<=max_mjd)
#         #        if 'exclude_bands' in kwargs:
#         #            for b in kwargs['exclude_bands']:
#         #                idx &= (d['band'] != b)
#         self.cache = d[idx]
#         self.f5_cache = self.f5[idx]
#         self.cache['mjd'] = np.floor(self.cache['mjd'])
#         self.cache_min_mjd = min_mjd
#         self.cache_max_mjd = max_mjd
        
#         # build a lighter structure, that contains just the visits 
#         logging.debug('cache fault -> updating visits [%7.0f,%7.0f]' % \
#                       (min_mjd, max_mjd))        
#         v = self.visits
#         idx = (v['mjd']>=min_mjd) & (v['mjd']<=max_mjd)
#         self.visit_cache = np.unique(v[idx])

#     def select_window(self, mjd, **kwargs):
#         """
#         Select the observation window, for a SN at a given redshift, that
#         peaks at a given mjd.
        
#         Args: 
#           mjd (float): SN peak 
#           z (float, optional): SN redshift [default 0.]
        
#         Returns:
#           r (ndarray of floats): 
#           u (ndarray of floats):
#         """
#         z = kwargs.get('z', 0.)
#         min_mjd, max_mjd = mjd + self.restframe_phase_range * (1.+z)
#         self._update_cache(min_mjd, max_mjd, **kwargs)
#         logging.debug('z: %5.2f mjd: [%f,%f]' % (z, min_mjd, max_mjd))        
        
#         # extract the cache contents that fit in the search window
#         idx = (self.cache['mjd']>min_mjd) & (self.cache['mjd']<max_mjd)
#         if 'band' in kwargs:
#             idx &= (self.cache['band'] == kwargs['band'])
#         r = self.cache[idx]
#         f5 = self.f5_cache[idx]
        
#         # extract the night cache contents that fit in the search window
#         idx = (self.visit_cache['mjd']>min_mjd) & (self.visit_cache['mjd']<max_mjd)
#         if 'band' in kwargs:    
#             idx &= (self.visit_cache['band'] == kwargs['band'])
#         u = self.visit_cache[idx]
        
#         return PixelLightcurveWindow(mjd, r, f5, u, etc=self.etc, unique_bandnames=self.unique_bandnames)

