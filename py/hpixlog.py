"""
"""

import logging 
logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)
import numpy as np
import healpy as hp


def extract_ddf_data(pxlobs, ddf_id=5, take_all_pixels=False):
    """Separate the DDF fields (healpixs) from the WFD fields (healpixs)
    """
    logging.info('separating the DDF fields (healpix) from the WFD/NGS/SGC')
    idx = pxlobs['proposal'] == ddf_id
    ddf_pixels = np.unique(pxlobs['pixel'][idx])    
    if not take_all_pixels:
        ddf_data = pxlobs[idx]
        wfd_data = pxlobs[~idx]
    else:
        logging.info('DDF pixels...')
        logging.info('%d DDF pixels identified' % len(ddf_pixels))
        idx = np.in1d(pxlobs['pixel'], ddf_pixels)
        ddf_data = pxlobs[idx]
        wfd_data = pxlobs[~idx]
    logging.info('WFD/NGS/SGC: %d obs, DDF: %d obs' % (len(wfd_data), len(ddf_data)))
    return wfd_data, ddf_data, ddf_pixels


def build_visits(data, exclude_bands=None):
    """Ligher index, containing pixel/band/night(mjd) information.
    
    Build a lighter structure, that contains just the information
    needed to count the number of individual nights during which a
    pixel was observed.
    """
    logging.info('build night/band array...')
    mjd = np.floor(data['mjd']).astype(int)
    visits = np.rec.fromarrays((data['pixel'].astype(int), mjd, data['band']),
                               names=['pixel', 'mjd', 'band'])
    logging.info('unique...')
    visits = np.unique(visits)
    logging.info('done...')
    
    idx = np.zeros(len(visits)).astype(bool)
    if exclude_bands is not None:
        for b in exclude_bands:
            idx |= (visits['band'] == b)
    logging.info('%d visits selected in bands: %r' % ((~idx).sum(), np.unique(visits[~idx]['band'])))
    return visits[~idx]


class PixelLightcurveWindow(object):
    """Manages the fraction around 
    
    """
    
    def __init__(self, peak_mjd, pixel_obs, f5, visits,
                 map_npix=hp.nside2npix(64),
                 min_restframe_phase_range=(-10, +30.),
                 max_restframe_phase_range=(-20, +45.), 
                 etc=None, unique_bandnames=None):
        self.peak_mjd = peak_mjd
        self.f5 = f5
        self.pixel_obs = pixel_obs
        self.visits = visits
        self.npix = map_npix
        self.min_restframe_phase_range = np.array(min_restframe_phase_range)
        self.max_restframe_phase_range = np.array(max_restframe_phase_range)
        self.etc = etc
        self.unique_bandnames = unique_bandnames

    def __len__(self):
        return len(self.visits)

    
    def compute(self, z, lcmodel):
        self.z = z
        self.snr, self.amplitude_snr = {}, {}
        for b in self.unique_bandnames:
            key = b + '-' + lcmodel.sntype
            self.amplitude_snr[key] = self.compute_amplitude_snr(z, b, lcmodel)
        key = lcmodel.sntype
        self.snr[key] = self.compute_snr(z, lcmodel)

    def first_last_visits(self, z=0.):
        """visits before and after a minimum restframe phase range. 
        """
        mjd_min_up,  mjd_max_low = self.peak_mjd + self.min_restframe_phase_range * (1.+z)
        mjd_min_low, mjd_max_up  = self.peak_mjd + self.max_restframe_phase_range * (1.+z)
        idx_min = (self.visits['mjd']>=mjd_min_low) & (self.visits['mjd']<=mjd_min_up)
        idx_max = (self.visits['mjd']>=mjd_max_low) & (self.visits['mjd']<=mjd_max_up)
        
        c_min = np.bincount(self.visits['pixel'][idx_min], minlength=self.npix).astype(float)
        c_max = np.bincount(self.visits['pixel'][idx_max], minlength=self.npix).astype(float)        
        return c_min, c_max

    def restframe_cadence(self, z, band=None):
        """observation cadence in the SN restframe
        """
        mjd_min, mjd_max = self.visits['mjd'].min(), self.visits['mjd'].max()
        dt = (mjd_max-mjd_min) / (1.+z)
        if band is None:
            c = np.bincount(self.visits['pixel'], minlength=self.npix).astype(float)
        else:
            idx = self.visits['band'] == band
            c = np.bincount(self.visits['pixel'][idx], minlength=self.npix).astype(float)
        c /= dt
        return c
        
    def compute_snr(self, z, lcmodel):
        """given a LC model and a redshift return the SNR of each point 
        """
        r = np.zeros(len(self.pixel_obs))
        d, f5 = self.pixel_obs, self.f5
        for b in self.unique_bandnames:
            idx = self.pixel_obs['bandname'] == b
            Li = lcmodel(d['mjd'][idx], self.peak_mjd, z, b)
            r[idx] = 5. * Li / f5[idx]
        
        return r
            
    def compute_amplitude_snr(self, z, bandname, lcmodel):
        """
        """
        c = np.zeros(self.npix)
        idx = self.pixel_obs['bandname'] == bandname
        d, f5 = self.pixel_obs[idx], self.f5[idx]
        Li = lcmodel(d['mjd'], self.peak_mjd, z, bandname)
        v = (Li/f5)**2
        snr = np.bincount(d['pixel'], weights=v, minlength=self.npix)
        return 5. * np.sqrt(snr)
        
    def compute_sigma_color(self, z, lcmodel):
        pass


class PixelObslog(object):
    """Iterable interface with a Pixel observation log 


    .. note : Since the typical pixel observation logs we deal with
          are large (10 years of survey => 2,4M exposures, nside=64 =>
          24,000 pixels in the southern sky, resulting in about 25M
          entries in the pixel log), this class manages a cache, in
          order to accelerate the requests.
    """
    
    def __init__(self, pxlobs, timestep=1, restframe_phase_range=(-25., 45.), 
                 max_redshift=1.2, exclude_bands=['u', 'y'], instrument_name='LSSTPG'):
        """
        """
        self.pixel_obs = self._clean_log(pxlobs,
                                         instrument_name=instrument_name,
                                         exclude_bands=exclude_bands)
        
        self.visits = build_visits(self.pixel_obs, exclude_bands=exclude_bands)
        self.mjd = np.unique(np.floor(pxlobs['mjd'])).astype(int)
        self.mjd = self.mjd[self.mjd>0] # why that ? 
        self.min_mjd, self.max_mjd = np.min(self.mjd), np.max(self.mjd)
        self.timestep = timestep
        self.cache_margin = np.array((-30., 100.))
        
        logging.info('loading ETC')
        self.etc = psf.find(instrument_name)
        self.f5 = self.etc.mag_to_flux(self.pixel_obs['m5'], self.pixel_obs['bandname'])
        
        self.restframe_phase_range = np.array(restframe_phase_range)
        self.sliding_window = np.array(restframe_phase_range) * (1. + max_redshift)
        self.max_redshift = max_redshift
        self.exclude_bands = exclude_bands
        logging.info('band names... ')
        self.unique_bandnames = np.unique(self.pixel_obs['bandname'])
        self.reinit()
        logging.info('ready to go')

    def _clean_log(self, data, instrument_name, exclude_bands):
        """
        """
        idx = np.zeros(len(data)).astype(bool)
        for b in exclude_bands:
            idx |= (data['band'] == b)
        logging.info('select: %d measurements out of %d' % ((~idx).sum(), len(idx)))
        d = data[~idx]
        
        suffix = instrument_name + '::'
        d['bandname'][:] = np.core.defchararray.add(suffix, d['band'])
        return d
        
    def reinit(self):
        """re-initialize the iterator to zero.
        """
        self.fig_odometer = 0
        self.cache_min_mjd = 0.
        self.cache_max_mjd = 0.
        self.cache = None
        self.current_mjd = self.min_mjd

    def __iter__(self):
        """iterator function. 
        """
        self.reinit()
        return self 

    def next(self):
        """the next function in the iterator 

        .. note : to be renamed __next__ when upgrading to python 3
        """
        dt = self.current_mjd + self.sliding_window
        self._update_cache(dt[0], dt[1], exclude_bands=self.exclude_bands)
        self.current_mjd += self.timestep
        r = self.select_window(self.current_mjd)
        if len(r) == 0 and (self.current_mjd > self.max_mjd):
            raise StopIteration
        return r
        
    # def _build_visits(self, exclude_bands=None):
    #     """Ligher index, containing pixel/band/night(mjd) information.
        
    #     Build a lighter structure, that contains just the information
    #     needed to count the number of individual nights during which a
    #     pixel was observed.
    #     """
    #     logging.info('build night/band array...')
    #     mjd = np.floor(self.pixel_obs['mjd']).astype(int)
    #     visits = np.rec.fromarrays((self.pixel_obs['pixel'].astype(int), mjd, self.pixel_obs['band']),
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
                
    def _update_cache(self, min_mjd, max_mjd, **kwargs):
        """Update the internal pixel log cache

        Since the pixel log is pretty large, we keep a smaller
        fraction of it, that is enough for ~50 iterations. This speeds
        up significantly the searches.
        
        """
        if self.cache is not None and \
           (max_mjd <= self.cache_max_mjd) and \
           (min_mjd >= self.cache_min_mjd):
            return 
        d = self.pixel_obs
        min_mjd = min_mjd + self.cache_margin[0]
        max_mjd = max_mjd + self.cache_margin[1]
        logging.debug('cache fault -> updating [%7.0f,%7.0f]' % \
                      (min_mjd, max_mjd))
        
        # select observations + cut on unused bands
        idx = (d['mjd']>=min_mjd) & (d['mjd']<=max_mjd)
        #        if 'exclude_bands' in kwargs:
        #            for b in kwargs['exclude_bands']:
        #                idx &= (d['band'] != b)
        self.cache = d[idx]
        self.f5_cache = self.f5[idx]
        self.cache['mjd'] = np.floor(self.cache['mjd'])
        self.cache_min_mjd = min_mjd
        self.cache_max_mjd = max_mjd
        
        # build a lighter structure, that contains just the visits 
        logging.debug('cache fault -> updating visits [%7.0f,%7.0f]' % \
                      (min_mjd, max_mjd))        
        v = self.visits
        idx = (v['mjd']>=min_mjd) & (v['mjd']<=max_mjd)
        self.visit_cache = np.unique(v[idx])

    def select_window(self, mjd, **kwargs):
        """
        Select the observation window, for a SN at a given redshift, that
        peaks at a given mjd.
        
        Args: 
          mjd (float): SN peak 
          z (float, optional): SN redshift [default 0.]
        
        Returns:
          r (ndarray of floats): 
          u (ndarray of floats):
        """
        z = kwargs.get('z', 0.)
        min_mjd, max_mjd = mjd + self.restframe_phase_range * (1.+z)
        self._update_cache(min_mjd, max_mjd, **kwargs)
        logging.debug('z: %5.2f mjd: [%f,%f]' % (z, min_mjd, max_mjd))        
        
        # extract the cache contents that fit in the search window
        idx = (self.cache['mjd']>min_mjd) & (self.cache['mjd']<max_mjd)
        if 'band' in kwargs:
            idx &= (self.cache['band'] == kwargs['band'])
        r = self.cache[idx]
        f5 = self.f5_cache[idx]
        
        # extract the night cache contents that fit in the search window
        idx = (self.visit_cache['mjd']>min_mjd) & (self.visit_cache['mjd']<max_mjd)
        if 'band' in kwargs:    
            idx &= (self.visit_cache['band'] == kwargs['band'])
        u = self.visit_cache[idx]
        
        return PixelLightcurveWindow(mjd, r, f5, u, etc=self.etc, unique_bandnames=self.unique_bandnames)
