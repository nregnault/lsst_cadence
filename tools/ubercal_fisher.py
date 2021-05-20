#!/usr/bin/env python 

import os
import os.path as op
import sys
import argparse

import logging
logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', 
                    level=logging.INFO)

import numpy as np
from scipy.sparse import coo_matrix, dia_matrix
from scipy.sparse import linalg
from sksparse import cholmod


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

import ubercal


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
    
    .. note: there is a degeneracy between the main zp's and 
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
    
    # if fix_grid_zp:
    #     idx = dp.gridobs == True
    #     tofix=np.unique(dp.expnum_index[idx])
    #     logging.info('%s zp parameters fixed: %r' %(len(tofix), tofix))
    #     model.params['zp'].fix(tofix, 0.)
    
    return model


def compute_block_hessian(model, y=None):
    class BlockHessian:
        pass
    ret = BlockHessian()
    
    logging.info('full jacobian')
    J = model.get_free_mat()

    logging.info('splicing jacobian')
    izp = model.params.indexof('zp') ; izp = izp[izp>=0]
    idzp = model.params.indexof('dzp') ; idzp = idzp[idzp>=0]
    Jp = J[:, list(izp) + list(idzp)]

    imag = model.params.indexof('mag') ; imag = imag[imag>=0]
    Jm = J[:,imag]
    
    logging.info('weight matrix')
    N = Jm.shape[0]
    s = 2.E-3
    wmap = (1./s)**2 * np.ones(N) 
    W = dia_matrix((wmap, [0]), shape=(N,N)) 
    
    logging.info('Hp')
    ret.Hp = Hp = Jp.T * W * Jp
    logging.info('Hm')
    ret.Hm = Hm = Jm.T * W * Jm
    logging.info('C')
    ret.C = C = Jp.T * W * Jm
    
    Hminv = Hm.tocoo()
    Hminv.data = 1. / Hminv.data
    ret.Hminv = Hminv
    
    logging.info('reduced upper zp matrix')
    ret.CHminv = CHminv = C * Hminv
    #    A = Hp - C * Hminv * C.T
    A = Hp - CHminv * C.T
    ret.A = A.tocsc()
    
    Bp, Bm = None, None
    if y is not None:
        Bp = Jp.T * W * y
        Bm = Jm.T * W * y
        
    ret.Bp = Bp
    ret.Bm = Bm

    return ret


def do_the_fit_by_blocks(bh, y):
    logging.info('cholesky: ')
    fact = cholmod.cholesky(bh.A)
    logging.info('calibration parameters')
    B = bh.Bp - bh.CHminv * bh.Bm
    p = fact(B)
    logging.info('star magnitudes')
    m = bh.Hminv * (bh.Bm - bh.C.T * p)
    logging.info('done')
    
    return p, m, fact


def save_fit_by_blocks_results(filename, p, m, dp, model, nside=1024):
    logging.info('fit results -> %s' % filename)
    N = hp.nside2npix(nside)
    mm = np.zeros(N)
    mm[dp.pixel_set] = m
    np.savez(filename, p=p, m=m, mm=mm)
    

def generate_fake_measurements(dp):
    y = np.zeros(len(dp.nt))
    dy = generate_gain_variations(dp)
    idx = dp.refstar == 1
    y[~idx] = -12. + dy[~idx]
    y[idx] = 19.
    return y


def getmat(model):
    """
    temporary overload of the getmat() method of linearmodel - with matrix dimensions explicitely defined.
    """
    N,n = len(model.rows), model.cols.max() + 1
    if model._valid is None:
        A = coo_matrix((model.vals, (model.rows, model.cols)), shape=(N,n)).tocsc()
    else:
        A = coo_matrix((model.vals[model._valid], (model.rows[model._valid], model.cols[model._valid])), shape=(N,n)).tocsc()
    return A

def get_free_mat(model):
    """
    temporary overload of the getmat() method of linearmodel - with matrix dimensions explicitely defined.
    """
    return getmat(model)[:, np.nonzero(model.params._free)[0]]


class FisherMatrix(object):
    """
    """
    def __init__(self, fn, refcell=94, nside=1024, perform_fit=False, refpixs=None):
        if type(fn) is str:
            d = [self._load_obslog(fn)]
        elif type(fn) is list:
            d = [self._load_obslog(f) for f in fn]
        if refpixs is not None:
            d += [self._add_refpixs(refpixs, d[0].dtype)]
        self.nside = nside
        self.obslog = np.hstack(d)
        #        self.obslog = self._remove_pixels_with_low_nbmeas(self.obslog, nmeas_min=5)
        self.npix = hp.nside2npix(nside)
        self.perform_fit = perform_fit
        self._flag_refcell_measurements(refcell)
        self._build_proxy()
        
    def _load_obslog(self, filename):
        logging.info('loading obslog: %s' % filename)
        f = np.load(filename)
        return f['l']
        
    def _flag_refcell_measurements(self, refcell):
        logging.info('identifying refcell measurements [%d]' % refcell)
        idx = self.obslog['cell'] == refcell
        self.obslog['refcell'][idx] = 1
        logging.info('%d measurements with that reference cell' % self.obslog['refcell'].sum())

    def _add_refpixs(self, refpixs, dtype):
        N = len(refpixs)
        z = np.zeros(N, dtype=dtype)
        z['pixel'] = refpixs
        z['refstar'] = 1
        return z

    def _remove_pixels_with_low_nbmeas(self, d, nmeas_min=5):
        logging.info('removing pixels with less observations than: %f' % nmeas_min)
        npix = hp.nside2npix(self.nside)
        c = np.bincount(d['pixel'], minlength=npix)
        pix = np.where(c<=nmeas_min)
        logging.info('')
        k = np.in1d(d['pixel'], pix)
        logging.info('%d pixels removed out of %d total measurements' % (k.sum(), len(d)))
        return d[k]
        
    def _build_proxy(self, focal_plane_map_validity=60):
        logging.info('building data proxy')
        dp = croaks.DataProxy(self.obslog, 
                              expnum='expnum', cell='cell', 
                              pixel='pixel', mjd='mjd', 
                              refstar='refstar',
                              refcell='refcell',
                              gridobs='gridobs')
        
        #
        # all of this may be moved to the model building functions
        # 
        # we need a value to index the focal plane validity periods
        mjd_min = dp.mjd[dp.mjd>0].min()
        s = (dp.mjd - mjd_min).astype(int) / focal_plane_map_validity
        s[s<0] = 0
        dp.add_field('month', s)
        
        # we need a unique value for the cell/months
        if dp.cell.max() < 1000:
            s = (dp.month * 1000 + dp.cell).astype(int)
        else:
            assert False
        dp.add_field('cellmonth', s)
        
        # indexes
        logging.info('building indexes')
        dp.make_index('pixel', intmap=1)
        dp.make_index('cell', intmap=1)
        dp.make_index('expnum', intmap=1)
        dp.make_index('cellmonth', intmap=1)

        # keep this temporarily (will be dropped later)
        self.dp = dp
        
        # the pixel set must survive the data proxy
        self.pixel_set = dp.pixel_set
        
        return dp

    def build_fisher_matrix(self, model, wmap=None):            
        logging.info('building the Fisher matrix')
        self.model = model
        self.params = self.model.params.copy()
        J = self.model.get_free_mat()
        N,n = J.shape
        
        if wmap is None:
            s = 2.E-3
            wmap = (1./s)**2 * np.ones(N)
        
        W = dia_matrix((wmap, [0]), shape=(N,N))
        JtW = J.T * W
        self.H = JtW * J
        #        self.H = J.T * W * J
        if self.perform_fit: 
            self.JtW = JtW
        logging.info('Fisher matrix: shape: %r, nnz=%d' % (self.H.shape, self.H.nnz))
        return self.H

    def cholesky(self, H):
        logging.info('cholesky factorization')
        self.fact = cholmod.cholesky(H)
        self.P = self.fact.P()
        return self.fact
                
    def low_mem_state(self):
        logging.info('dropping now useless data structures')
        del self.obslog
        del self.model 
        
    def selinv(self):
        """
        compute the diagonal of the inverse Fisher matrix 
        (we use selinv and we recycle the permulation given by cholmod)
        """
        H, P = self.H, self.P
        logging.info('diagonal of the inverse Fisher matrix')
        logging.info('shapes: H: %r - nnz=%d' % (H.shape,H.nnz))
        #        logging.info('shapes: P: %r ' % (P.shape,))
        
        r = saunerie.selinv.selinv(H, P=P)
        md = np.zeros(hp.nside2npix(nside))
        md[:] = hp.UNSEEN
        p = model.params.copy()
        p.free = r
        md[self.pixel_set] = np.sqrt(p['mag'].full)
    
        return md
        
    def survey_realizations(self, n_realizations=1):
        ret = []
        D = self.fact.D()
        p = self.params.copy()
        
        for i_realization in xrange(n_realizations):
            logging.info('survey realizations: %d : cell dmags' % i_realization)
            z = np.random.normal(size=len(p.free))
            zz = z / np.sqrt(D)
            zz = self.fact.solve_Lt(zz)
            zz = self.fact.apply_Pt(zz)
            m = np.zeros(hp.nside2npix(self.nside))
            m[:] = hp.UNSEEN
            p.free = zz
            m[self.pixel_set] = p['mag'].full
            ret.append(m)
        return ret
        
    def fit(self, y):
        logging.info('survey fit:')
        p = self.params.copy()
        m = np.zeros(self.npix)
        m[:] = hp.UNSEEN
        
        B = self.JtW * y
        p.free = self.fact(B)
        m[self.pixel_set] = p['mag'].full
        return p,m
        

def generate_gain_variations(dp, slopes=None, dt=4/24.):
    if slopes is None:
        logging.info('generating slopes:')
        cellid = np.unique(dp.cell)
        slopes = np.random.uniform(0., 0.005, size=len(cellid))
    # dates within the night 
    logging.info('generation delta meas. due to gain variations')
    dt = 4. / 24.
    m = np.floor(dp.mjd - dt) 
    dt_vs_midnight = dp.mjd % m - dt - 0.5
    dy = slopes[dp.cell] * dt_vs_midnight
    return dy 


def compute_cl(m, nside=1024, nest=True, lmax=1500):
    logging.info('compute_cl -> lmax=%d' % lmax)
    if nest:
        npix = hp.nside2npix(nside)
        mm = np.zeros(npix)
        i = np.arange(npix)
        mm[:] = m[hp.ring2nest(nside,i)]
    else:
        mm = m
    c = hp.anafast(mm, lmax=lmax)
    return c


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""compute and factorize the fisher matrix""")
    parser.add_argument('-O', '--output_dir',
                        default=None,
                        help='output directory')
    parser.add_argument('--nside', dest='nside', default=1024, type=int,
                        help='HEALPIX pixellization')
    # parser.add_argument('--ncellsx', default=1, type=int, dest='nx',
    #                     help='number of cells per ccd in the x-direction')
    # parser.add_argument('--ncellsy', default=1, type=int, dest='ny',
    #                     help='number of cells per ccd in the y-direction')
    parser.add_argument('--dump_inv_fisher', default=None, type=str,
                        help='dump the fisher matrix (in coo_matrix format) in specified dir')
    parser.add_argument('--selinv-for-inversion', default=False,
                        action='store_true',
                        help='use selinv to get the diagonal errors (takes time, default is MC inversion)')
    parser.add_argument('--dump_maps', default=False, 
                        action='store_true',
                        help='dump the survey mag error realizations (takes disk space)')
    parser.add_argument('--dump_cls', default=False, 
                        action='store_true',
                        help="dump the map Cl's")
    parser.add_argument('--realizations', dest='realizations',
                        default=100, type=int, 
                        help='number of realizations')
    parser.add_argument('--lmax', dest='lmax',
                        default=1500, type=int, 
                        help='lmax for cl')
    parser.add_argument('--plot_dir', default=None, type=str,
                        help='prepare and dump control plots')
    parser.add_argument('--fit', default=0, 
                        help='perform fits')
    parser.add_argument('--refpixs', type=str, default='equatorial',
                        help='where to place refpixels [equatorial|ddf_4|ddf_all|random]')
    parser.add_argument('obs', type=str, nargs='+',
                        help='observation file(s)')
    
    args = parser.parse_args()
    print args    

    # reference pixels 
    if args.refpixs == 'equatorial':
        refpixs = ubercal.get_flux_standards(args.nside, equatorial=True)
        logging.info('equatorial refs: refpixs=%r' % refpixs)
    elif args.refpixs == 'ddf_4':
        refpixs = ubercal.get_flux_standards(args.nside, ddf_4=True)
        logging.info('DDF refs: refpixs=%r' % refpixs)
    elif args.refpixs == 'ddf_all':
        refpixs = ubercal.get_flux_standards(args.nside, ddf_all=True)
        logging.info('DDF refs: refpixs=%r' % refpixs)
    elif random:
        refpixs = ubercal.get_flux_standards(args.nside, random=10)
        logging.info('random refs: refpixs=%r' % refpixs)
    else:
        refpixs = ubercal.get_flux_standards(args.nside)
        logging.info('default refs: refpixs=%r' % refpixs)

    
    # compute the fisher matrix
    #    dp, model, H, fact = fisher_matrix(args.obs)
    fm = FisherMatrix(args.obs, refcell=94, perform_fit=True, refpixs=refpixs)
    
