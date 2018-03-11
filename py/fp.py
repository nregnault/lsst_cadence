"""
"""


import numpy as np
from scipy.sparse import coo_matrix, dia_matrix
from scipy.sparse import linalg
from sksparse import cholmod

import logging 
logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', 
                    level=logging.INFO)

import matplotlib
# matplotlib.use('Agg')

from matplotlib.patches import Rectangle, Polygon
from matplotlib.collections import PatchCollection
import pylab as pl
pl.interactive(0)

import healpy as hp
from mx.DateTime import DateTimeFrom

import croaks
import saunerie
from saunerie.fitparameters import FitParameters 
import saunerie.selinv as selinv
import saunerie.linearmodels as lm



RADIAN_PER_DEGREE = np.pi / 180.
DEGREE_PER_RADIAN = 180. / np.pi


            
            



def main_map(fn='all_fields.npy'):
    nside = 1024
    m = np.zeros(hp.nside2npix(nside))    

    # opsim log 
    log = np.load(fn)
    mjd_min = DateTimeFrom('2022-01-01').mjd
    mjd_max = DateTimeFrom('2023-01-01').mjd
    idx = (log['band'] == 'LSSTPG::i')
    idx &= (log['mjd'] >= mjd_min)
    idx &= (log['mjd'] <= mjd_max)
    log = log[idx]
    coords = np.vstack((log['Ra'] * DEGREE_PER_RADIAN, log['Dec'] * DEGREE_PER_RADIAN)).T
    print coords.shape 
        
    f = FocalPlane()
    p = f.pixellize(1,1)    
    
    for i,(r,d) in enumerate(coords):
        if i%100 == 0:
            print i, r, d
        cells = p.copy()
        f.to_uv(cells, (r,d))
        
        for c in cells:
            poly = np.vstack([c['p0'], c['p1'], c['p2'], c['p3']])
            ipix = hp.query_polygon(nside, poly, nest=True)
            m[ipix] += 1.
            
    return m, p


def query(cells, nside=1024):
    for c in cells:
        poly = np.vstack([c['p0'], c['p1'], c['p2'], c['p3']])
        ipix = hp.query_polygon(nside, poly, nest=True)


def build_proxy(l):
    dp = croaks.DataProxy(l, expnum='expnum', cell='cell', 
                          pixel='pixel', mjd='mjd', 
                          refstar='refstar',
                          refcell='refcell',
                          gridobs='gridobs')
    
    mjd_min = dp.mjd[dp.mjd>0].min()
    s = (dp.mjd - mjd_min).astype(int) / 60
    s[s<0] = 0
    dp.add_field('month', s)
    
    if dp.cell.max() < 1000:
        s = (dp.month * 1000 + dp.cell).astype(int)
    else:
        assert False
    dp.add_field('cellmonth', s)
    dp.make_index('pixel', intmap=1)
    dp.make_index('cell', intmap=1)
    dp.make_index('expnum', intmap=1)
    dp.make_index('cellmonth', intmap=1)
    
    return dp


def simple_ubercal_model(dp):
    """
    the simplest possible ubercal model: 
    """
    mag = lm.indic(dp.pixel_index, name='mags')
    zpexp = lm.indic(dp.expnum_index, name='zps', valid=dp.refstar==0)
    return zpexp + mag


def ubercal_model_with_photflat(dp, fix_grid_zp=True):
    """
    Ubercal model with one photometric flat on
    
    .. note : there is a degeneracy between the main zp's and 
              the photometric flat dzp's and one needs to fix
              at least one of them per epoch.
    
    """
    mag = lm.indic(dp.pixel_index, name='mag')
    zp = lm.indic(dp.expnum_index, name='zp', valid=dp.refstar==0)
    dzp = lm.indic(dp.cellmonth_index, name='dzp', valid=dp.refstar==0)
    model = mag + zp + dzp
    
    idx = dp.refcell == 1
    tofix = np.unique(dp.cellmonth_index[idx])
    logging.info('%s dzp parameters fixed: %r' %(len(tofix), tofix))
    model.params['dzp'].fix(tofix, 0.)
    #    model.params['dzp'].fix()
    
    if fix_grid_zp:
        idx = dp.gridobs == True
        tofix=np.unique(dp.expnum_index[idx])
        logging.info('%s zp parameters fixed: %r' %(len(tofix), tofix))
        model.params['zp'].fix(tofix, 0.)
    
    return model


class UbercalModel(object):
    """
    I no longer use that. 
    """
    def __init__(self, data):
        self.data = data
        # identify the N (pseudo-)stars that have measurements
        # and remap their indices to [0;N-1]
        # self.pixels contains the (unique) healpix index of the pixels 
        # self.istars contains the remapped values of these indices
        logging.info('[UbercalModel] identify measured stars')
        self.pixels, self.istar = np.unique(data['pixel'], return_inverse=1)
        self.nstars = len(self.pixels)
        
        # exposures 
        logging.info('[UbercalModel] identify exposures')
        exp, self.iexp = np.unique(data['expnum'], return_inverse=1)
        self.nexp = len(exp)
        # internal parameter vector 
        self.pars = self.init_pars()
        
    def init_pars(self):
        p = FitParameters([('zps', self.nexp), ('mags', self.nstars)])
        return p

    def __call__(self, p=None, jac=False):
        
        if p is not None:
            self.pars.free = p
        
        m = self.pars['mags'].full[self.istar]
        zp = self.pars['zps'].full[self.iexp]
        val = m + zp
        
        N = len(self.data)
        n = len(self.pars.free)
        nz = len(self.pars['zps'].full)        
        
        if jac:
            i_zp, v_zp = np.arange(N), np.ones(N, dtype='f8')
            j_zp_free = self.pars['zps'].indexof(self.iexp)
            idx = j_zp_free >= 0
            i_zp, j_zp, v_zp = i_zp[idx], j_zp_free[idx], v_zp[idx]
            
            i_mags, v_mags = np.arange(N), np.ones(N, dtype='f8')
            j_mags_free = self.pars['mags'].indexof(self.istar)
            idx = j_mags_free >= 0
            i_mags, j_mags, v_mags = i_mags[idx], j_mags_free[idx], v_mags[idx]
            
            i = np.hstack((i_zp, i_mags))
            j = np.hstack((j_zp, j_mags))
            v = np.hstack((v_zp, v_mags))
            J = coo_matrix((v, (i,j)), shape=(N,n))
            return val, J

        return val

        

def main_observe(nside=1024, refpixs=None, fn=None, dither=False, mjd_range=None, band='LSSTPG::i'):
    """
    Build the observation NTuple
    (each measurement associates a focal plane cell to a sky pixel)
    """
    m = np.zeros(hp.nside2npix(nside))

    # opsim log 
    log = np.load(fn)
    if mjd_range is None:
        mjd_min = DateTimeFrom('2023-01-01').mjd
        mjd_max = DateTimeFrom('2024-01-01').mjd
    else:
        mjd_min, mjd_max = mjd_range
    idx = (log['band'] == band)
    idx &= (log['mjd'] >= mjd_min)
    idx &= (log['mjd'] <= mjd_max)
    log = log[idx]
    iexp = np.arange(len(log))
    n_exp = len(log)                                                                                                                                                                          
    if dither:
        d_ra = np.random.uniform(-1. * RADIAN_PER_DEGREE, 1. * RADIAN_PER_DEGREE, n_exp) 
        d_dec = np.random.uniform(1. * RADIAN_PER_DEGREE, 1. * RADIAN_PER_DEGREE, n_exp)
    else:
        d_ra = np.zeros(n_exp)
        d_dec = np.zeros(n_exp)
        
    coords = np.vstack((log['Ra'] * DEGREE_PER_RADIAN, 
                        log['Dec'] * DEGREE_PER_RADIAN, 
                        log['mjd'],
                        iexp,
                        log['gridobs'],
                        log['rotation'])).T
    dtype=np.dtype([('expnum', 'i4'), ('cell', 'i2'), 
                    # ('ichip', 'i4'), 
                    # ('iraft', 'i4'), 
                    ('pixel', 'i4'), ('mjd', 'f8'),
                    ('refstar', 'i2'),
                    ('refcell', 'i2'),
                    ('gridobs', 'i2')])
    
    f = FocalPlane()
    p = f.pixellize(1,1)    
    l = []
    
    for i,(r,d,mjd,iexp,gridobs,rotation) in enumerate(coords):
        if i%100 == 0:
            logging.info('exposure: % 6d  [RA=% 6.6f DEC=% +6.5f]' % (i,r,d))
        cells = p.copy()
        if np.abs(rotation)>0.:
            cells = f.rotate(cells, theta=rotation * 180. / np.pi)
        f.to_uv(cells, (r,d))
        
        for c in cells:
            poly = np.vstack([c['p0'], c['p1'], c['p2'], c['p3']])
            try:
                ipix = hp.query_polygon(nside, poly, nest=True, inclusive=False)
            except:
                continue
            m[ipix] += 1
            N = len(ipix)
            z = np.zeros(N, dtype=dtype)
            z['expnum'] = iexp
            z['cell'] = c['icell']
            z['pixel'] = ipix
            z['mjd'] = mjd
            z['gridobs'] = gridobs
            l.append(z)

    # reference objects (well, here, we have reference pixels)
    if refpixs is not None:
        N = len(refpixs)
        z = np.zeros(N, dtype=dtype)
        z['pixel'] = refpixs
        z['refstar'] = 1
        l.append(z)
        
    return m, np.hstack(l)


def radec2pix(ra_deg, dec_deg, nside=1024):
    phi   = ra_deg * RADIAN_PER_DEGREE
    theta = 0.5*np.pi - dec_deg * RADIAN_PER_DEGREE
    return hp.ang2pix(nside, theta, phi, nest=True)


def get_flux_standards(nside, random=0):
    if random <= 0.:
        ra = np.array([0., 12., 32., 55.])
        dec = np.array([0., -22., -12., -45.])
    else:
        ra = np.random.uniform(0., 350., random)
        dec = np.random.uniform(-55, -5., random)
    return radec2pix(ra, dec, nside=nside)


def plot_flux_standards(ra, dec, nside=256):
    #    ra = np.array([0., 12., 32., 55.])
    #    dec = np.array([0., -22., -12., -45.])
    pix = radec2pix(ra, dec, nside=nside)
    m = np.zeros(hp.nside2npix(nside))
    m[pix] = 1.E6
    hp.mollview(m, nest=True, xsize=5000)
    

def main_fisher(fn='tst.obs.npz', refcell=94):
    logging.info('loading data from: %s' % fn)
    t = np.load(fn)
    m,l = t['m'], t['l']
    
    logging.info('identifying refcell measurements [%d]' % refcell)
    idx = l['cell'] == refcell
    l['refcell'][idx] = 1
    logging.info('%d measurements on that reference cell' % l['refcell'].sum())
    
    logging.info('data proxy')
    dp = build_proxy(l)
    
    logging.info('ubercal model')
    #    model = simple_ubercal_model(dp)
    model = ubercal_model_with_photflat(dp)

    logging.info('building the Fisher matrix')
    J = model.get_free_mat()
    N,n = J.shape
    W = dia_matrix((1.E5 * np.ones(N), [0]), shape=(N,N))
    H = J.T * W * J
    
    logging.info('cholesky factorization')
    fact = cholmod.cholesky(H)
    
    return dp, model, H, fact
    

def main_selinv(dp, model, H, fact, nside=1024, plot=True):
    # compute and plot the diagonal of the inverse Fisher matrix
    logging.info('diagonal of the inverse Fisher matrix')
    r = saunerie.selinv.selinv(H, fact.P())
    md = np.zeros(hp.nside2npix(nside))
    md[:] = hp.UNSEEN
    p = model.params.copy()
    p.free = r
    md[dp.pixel_set] = np.sqrt(p['mag'].full)
    
    if plot is True:
        hp.mollview(md, nest=True)
        
    return md

def main_realizations(dp, model, H, fact, nside=1024, n_realizations=1, plot=True, dump_plot=False, prefix='realization'):
    if dump_plot == True:
        plot = True
    # compute and plot the specified number of realizations
    # of the fitted mags
    ret = []
    D = fact.D()
    for i_realization in xrange(n_realizations):
        logging.info('survey realizations: %d : cell dmags' % i_realization)
        z = np.random.normal(size=len(model.params.free))
        zz = fact.solve_Lt(z)
        zz = zz / np.sqrt(D)
        zz = fact.apply_Pt(zz)
        m = np.zeros(hp.nside2npix(nside))
        m[:] = hp.UNSEEN
        p = model.params.copy()
        p.free = zz
        m[dp.pixel_set] = p['mag'].full
        if plot is True:
            hp.mollview(m, nest=True, min=-0.01, max=0.01)
            if dump_plot:
                pl.gcf().savefig(prefix + "_%03d.png" % i_realization)
        ret.append(m)
    return ret


def stack_cl(dp, model, H, fact, nside=1024, n_realizations=500):
    cl = []
    npix = hp.nside2npix(nside)
    mm = np.zeros(npix)
    i = np.arange(npix)
    for i_realization in xrange(n_realizations):
        ret = main_realizations(dp, model, H, fact, nside=nside, n_realizations=1, plot=False)
        m = ret[0]
        mm[:] = m[hp.ring2nest(1024,i)]
        c = hp.anafast(mm, lmax=1500)
        cl.append(c)
        del ret
    return cl


def observe_all(pkg=0):
    if pkg == 0:
        m,l = main_observe(refpixs=get_flux_standards(1024), dither=False, fn='all_fields_wfd_dd_minion1626_grid_25_rotation.npy')
        np.savez('obs_minion1626_grid_25_rotation.npz', l=l, m=m)
        hp.mollview(m, nest=1)
        pl.gcf().savefig('nobs_minion1626_grid_25_rotation.png')
    elif pkg == 1:
        m,l = main_observe(refpixs=get_flux_standards(1024), dither=False, fn='all_fields_wfd_dd_minion1626_grid_25_norotation.npy')
        np.savez('obs_minion1626_grid_25_norotation.npz', l=l, m=m)
        hp.mollview(m, nest=1)
        pl.gcf().savefig('nobs_minion1626_grid_25_norotation.png')
    elif pkg == 2:
        m,l = main_observe(refpixs=get_flux_standards(1024), dither=False, fn='all_fields_wfd_dd_minion1626_nogrids.npy')
        np.savez('obs_minion1626_nogrids.npz', l=l, m=m)
        hp.mollview(m, nest=1)
        pl.gcf().savefig('nobs_minion1626_nogrids.png')
    elif pkg == 3:
        m,l = main_observe(refpixs=get_flux_standards(1024), dither=False, fn='all_fields_wfd_dd_baseline10yr_grid_25_norotation.npy')
        np.savez('obs_baseline_grid_25_norotation.npz', l=l, m=m)
        hp.mollview(m, nest=1)
        pl.gcf().savefig('nobs_baseline_grid_25_norotation.png')
    elif pkg == 4:
        m,l = main_observe(refpixs=get_flux_standards(1024), dither=False, fn='all_fields_wfd_dd_baseline10yr_nogrids.npy')
        np.savez('obs_baseline_nogrids.npz', l=l, m=m)
        hp.mollview(m, nest=1)
        pl.gcf().savefig('nobs_baseline_nogrids.png')
        
#all_fields_wfd_dd_baseline10yr_nogrids.npy
#all_fields_wfd_dd_baseline10yr_grid_25_norotation.npy

##################################################################################
def main_fisher_diag(l, neigs=0):
    
    logging.info('main_fisher - new UbercalModel')
    model = UbercalModel(l)
    logging.info('building fisher matrix')
    model.pars['zps'].fix(0, 0.0)
    v,J = model(jac=True)
    #    J.data /= 0.01**2
    N,n = J.shape
    H = J.T * J 
    
    logging.info('cholesky factorization')
    fact = cholmod.cholesky(H)
    
    logging.info('cholesky inverse')
    def mv(x):
        return fact.solve_A(x)
    OPinv = OPinv = linalg.LinearOperator(H.shape, mv)
    
    e,v = None, None
    if neigs > 0:
        logging.info('extracting smallest eigenvals/eigenvects')
        #    X = np.random.random((n,3))
        X = np.random.random(neigs * 5288154)
        X = X = X.reshape((-1,neigs))
        tol = 1.E-5
        e,v = linalg.lobpcg(H,X, largest=False, verbosityLevel=2, maxiter=100, M=OPinv)
    
    return H, e, v, model, fact


def draw_map(m, cmap=pl.cm.jet, filename=None):
    mm = m.copy()
    mm[mm<=0] = hp.UNSEEN
    cmap.set_under('w')
    hp.mollview(mm, nest=True, cmap=cmap)
    hp.graticule()
    if filename is not None:
        fig = pl.gcf()
        fig.savefig(filename + '.png', bbox_inches='tight')


def draw_eigenvect(i, e, v, model, cmap=pl.cm.jet, filename=None, sig=None):
    if sig is not None:
        scale = np.sqrt(sig**2 / e[i])
    else:
        scale = 1.
    print scale 
    m = np.zeros(hp.nside2npix(1024))
    m[:] = hp.UNSEEN
    vv = v[:,i][model.pars['mags'].indexof()]
    cmap.set_under('w')
    m[model.pixels] = scale * vv
    hp.mollview(m, nest=True, cmap=cmap)
    hp.graticule()
    if filename is not None:
        fig = pl.gcf()
        fig.savefig(filename + '.png', bbox_inches='tight')
    return m
    

    
    
def extreme_eigenvecs(model, dp, H, fact, nside=1024):
    N = H.shape[0]    
    logging.info('eigenector corresponding to the smallest eigen value of C')
    X = np.random.random(size=N)
    for i in xrange(2):
        X /= np.sqrt((X**2).sum())
        X = H * X
    
    m = np.zeros(hp.nside2npix(nside))
    m[:] = hp.UNSEEN
    X /= np.sqrt((X**2).sum())
    model.params.free = X
    m[dp.pixel_set] = model.params['mag'].full
    hp.mollview(m, nest=True)
    
    logging.info('eigenector corresponding to the largest eigen value of C')
    X = np.random.random(size=N)
    for i in xrange(5):
        X /= np.sqrt((X**2).sum())
        Y = fact.solve_LDLt(X)
        X = Y.copy()
    
    m = np.zeros(hp.nside2npix(nside))
    m[:] = hp.UNSEEN
    X /= np.sqrt((X**2).sum())
    model.params.free = fact.apply_Pt(X)
    m[dp.pixel_set] = model.params['mag'].full # X[model.pars['mags'].indexof()]
    hp.mollview(m, nest=True)
    

# def main_fisher(dp, model, nside=1024):
#     logging.info('main_fisher - new UbercalModel')
#     #    model = UbercalModel(l)
#     logging.info('building the Fisher matrix')
#     #    model.pars['zps'].fix(0, 0.0)
#     #    for i in xrange(10):
#     #        model.pars['mags'].fix(42 + 1000*i, 0.0)
#     v,J = model(jac=True)
    
#     #    J.data /= 0.01**2
#     N,n = J.shape

#     # nothing very realistic here. 
#     # just to rescale the uncertainties
#     W = dia_matrix((1.E4 * np.ones(N), [0]), shape=(N,N))
#     H = J.T * W * J 
    
#     logging.info('cholesky factorization')
#     fact = cholmod.cholesky(H)

#     return model, H, fact


# def main_old():

#     ra, dec = -45., 12.
#     nside = 1024
#     m = np.zeros(hp.nside2npix(nside))
#     f = FocalPlane()    
#     p = f.pixellize(1,1)
    
#     cells = p.copy()
#     f.to_uv(cells, (ra,dec))
    
#     for c in cells:
#         poly = np.vstack([c['p0'], c['p1'], c['p2'], c['p3']])
#         ipix = hp.querypolygon(nside, poly)
#         m[ipix] += c['ichip']
    
#     return m, cells

    
# def main_random(N=100):
    
#     nside = 1024
#     m = np.zeros(hp.nside2npix(nside))    
    
#     f = FocalPlane()
#     p = f.pixellize(1,1)    
#     coords = np.vstack((ra,dec)).T
    
#     for i,(r,d) in enumerate(coords):
#         if i%100 == 0:
#             print i, r, d
#         cells = p.copy()
#         f.to_uv(cells, (r,d))
        
#         for c in cells:
#             poly = np.vstack([c['p0'], c['p1'], c['p2'], c['p3']])
#             ipix = hp.query_polygon(nside, poly, nest=True)
#             m[ipix] += 1. # c['icell']
            
#     return m, p


# def draw_map(cells, lon_0=0., lat_0=-70):
#     fig = pl.figure(figsize=(8,8))
#     #    m = Basemap(projection='hammer', lon_0=0.)
#     m = Basemap(projection='gnom', 
#                 lon_0=lon_0, lat_0=lat_0,
#                 llcrnrlon=lon_0-30., llcrnrlat=lat_0-30.,
#                 urcrnrlon=lon_0+30., urcrnrlat=lat_0+30.,)
#     m.drawmeridians(np.arange(lon_0, lon_0+360., 10.))
#     m.drawparallels(np.arange(lat_0-90., lat_0+90., 10.))
#     x,y = m(cells['p0'][:,0] * DEGREE_PER_RADIAN, cells['p0'][:,1] * DEGREE_PER_RADIAN)
#     m.plot(x,y, 'r.')


