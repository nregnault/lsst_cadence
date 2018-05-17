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


class FisherMatrix(object):
    """
    """
    def __init__(self, filename, refcell=94, nside=1024):
        self.f = np.load(filename)
        self.nside = nside
        self.obslog = self._load_obslog(self.f)
        self.flag_refcell_measurements(refcell)
        self._build_proxy()
        
    def _load_obslog(self):
        logging.info('loading obslog: %s' % self.filename)
        return self.f['l']
        
    def _flag_refcell_measurements(self, refcell):
        logging.info('identifying refcell measurements [%d]' % refcell)
        idx = self.obslog['cell'] == refcell
        self.obslog['refcell'][idx] = 1
        logging.info('%d measurements with that reference cell' % self.obslog['refcell'].sum())
        
    def _build_proxy(self, focal_plane_map_validity=30):
        logging.info('building data proxy')
        dp = croaks.DataProxy(l, expnum='expnum', cell='cell', 
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
            wmap = 1.E6 * np.ones(N)
        
        W = dia_matrix((wmap, [0]), shape=(N,N))
        self.H = J.T * W * J
        return self.H

    def cholesky(self, H):
        logging.info('cholesky factorization')
        self.fact = cholmod.cholesky(H)
        self.P = self.fact.P()
        return self.fact
        
    def low_mem_state(self):
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
        logging.info('shapes: P: %r ' % (P.shape,))
        
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
            zz = fact.solve_Lt(zz)
            zz = fact.apply_Pt(zz)
            m = np.zeros(hp.nside2npix(self.nside))
            m[:] = hp.UNSEEN
            p.free = zz
            m[self.pixel_set] = p['mag'].full
            ret.append(m)
        return ret
        


def compute_cl(m, nside=1024, nest=True, lmax=1500):
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
    parser.add_argument('obs', type=str, default=None,
                        help='observation file')
    
    args = parser.parse_args()
    print args    
    
    # compute the fisher matrix
    #    dp, model, H, fact = fisher_matrix(args.obs)
    fm = FisherMatrix(args.obs, refcell=94)
    
    logging.info('ubercal model')
    #    model = simple_ubercal_model(fm.dp)
    model = ubercal_model_with_photflat(fm.dp)
    fm.build_fisher_matrix(model)
    fm.low_mem_state()    
    fm.cholesky(fm.H)
    
    if args.dump_inv_fisher:
        fn = args.dump_inv_fisher + os.sep + 'H.npz'
        HH = fm.H.tocoo()
        np.savez(fn, v=HH.data, i=HH.row, j=HH.col, permutation=fact.P())
        
    # get the diagonal errors using selinv
    # takes about as much time as the cholesky factorization itself
    if args.selinv_for_inversion:
        md = fm.selinv()
        if args.output_dir is not None:
            np.save(args.output_dir + os.sep + 'md.npy', md)
        if args.plot_dir is not None:
            hp.mollview(md, nest=True, min=0., max=0.001)
            pl.gcf().savefig(args.plot_dir + os.sep + 'md.png', bbox_inches='tight')

    if args.realizations > 0:
        # by default, we generate a handful of realizations, and perform
        sr = fm.survey_realizations(n_realizations=args.realizations)
        sr = np.vstack(sr)
        cl = [compute_cl(r, nside=args.nside, nest=True, lmax=args.lmax) for r in ret]
        cl = np.vstack(cl)

        # dumps the maps and the cl's (optional)
        if args.output_dir is not None:
            if args.dump_maps:
                np.save(args.output_dir + os.sep + 'survey_realizations.npy', sr)
            if args.dump_cls:
                np.save(args.output_dir + os.sep + 'cl.npy', cl)
                
        if args.plot_dir is not None:
            for i,r in enumerate(sr):
                hp.mollview(r, nest=True, min=-0.001, max=0.001)
                pl.gcf().savefig(args.plot_dir + os.sep + 'r_%05d.png' % i, bbox_inches='tight')
                
            fig = pl.figure()
            l = np.arange(1, args.lmax+2, 1.)
            pl.plot(l*(l+1)*cl[0], 'k,-')
            fig.savefig(args.plot_dir + os.sep + 'cl_0.png', bbox_inches='tight')
            fig = pl.figure()
            for c in cl:
                pl.plot(l*(l+1)*c.mean(axis=0), 'k,-')
                fig.savefig(args.plot_dir + os.sep + 'cl_all.png', bbox_inches='tight')
        
    



# def build_proxy(l):
#     dp = croaks.DataProxy(l, expnum='expnum', cell='cell', 
#                           pixel='pixel', mjd='mjd', 
#                           refstar='refstar',
#                           refcell='refcell',
#                           gridobs='gridobs')
    
#     mjd_min = dp.mjd[dp.mjd>0].min()
#     s = (dp.mjd - mjd_min).astype(int) / 15
#     s[s<0] = 0
#     dp.add_field('month', s)
    
#     if dp.cell.max() < 1000:
#         s = (dp.month * 1000 + dp.cell).astype(int)
#     else:
#         assert False
#     dp.add_field('cellmonth', s)
#     dp.make_index('pixel', intmap=1)
#     dp.make_index('cell', intmap=1)
#     dp.make_index('expnum', intmap=1)
#     dp.make_index('cellmonth', intmap=1)
    
#     return dp
        

# parser.add_argument('--mjd_min',
#                     dest='mjd_min', default=59580., type=float,
#                     help='observation period start')
# parser.add_argument('--mjd_max',
#                     dest='mjd_max', default=59945., type=float,
#                     help='observation period end')
# parser.add_argument('--band',
#                     dest='band', default='z', type=str,
#                     help='LSST passband [ugrizy]')



# def fisher_matrix(fn='tst.obs.npz', refcell=94):
#     """
#     Compute the fisher matrix from an observation log.
#     """
#     logging.info('loading data from: %s' % fn)
#     t = np.load(fn)
#     m,l = t['m'], t['l']
    
#     logging.info('identifying refcell measurements [%d]' % refcell)
#     idx = l['cell'] == refcell
#     l['refcell'][idx] = 1
#     logging.info('%d measurements on that reference cell' % l['refcell'].sum())
    
#     logging.info('data proxy')
#     dp = build_proxy(l)
    
#     logging.info('ubercal model')
#     #    model = simple_ubercal_model(dp)
#     model = ubercal_model_with_photflat(dp)

#     logging.info('building the Fisher matrix')
#     J = model.get_free_mat()
#     N,n = J.shape
#     W = dia_matrix((1.E6 * np.ones(N), [0]), shape=(N,N))
#     H = J.T * W * J
    
#     logging.info('cholesky factorization')
#     fact = cholmod.cholesky(H)
    
#     return dp, model, H, fact
    


# def fisher_selinv(dp, model, H, fact, nside=1024):
#     """
#     compute the diagonal of the inverse Fisher matrix 
#     (we use selinv and we recycle the permulation given by cholmod)
#     """
#     logging.info('diagonal of the inverse Fisher matrix')
#     logging.info('shapes: H: %r - nnz=%d' % (H.shape,H.nnz))
#     logging.info('shapes: P: %r ' % (fact.P().shape,))
#     logging.info('debug: P: (%E,%E)' % (fact.P().min(), fact.P().max()))
#     r = saunerie.selinv.selinv(H, P=fact.P())
#     #    r = saunerie.selinv.selinv(H)
#     md = np.zeros(hp.nside2npix(nside))
#     md[:] = hp.UNSEEN
#     p = model.params.copy()
#     p.free = r
#     md[dp.pixel_set] = np.sqrt(p['mag'].full)
    
#     return md


#    # genrate a handful of realizations 
#    ret = generate_realizations(dp, model, H, fact, 
#                                nside=args.nside, 
#                                n_realizations=args.realizations)




# def generate_realizations_from_fisher_matrix(dp, model, H, fact, nside=1024, n_realizations=1):
#     """
#     Generate one realization of the survey errors 
#     """
#     ret = []
#     D = fact.D()
#     for i_realization in xrange(n_realizations):
#         logging.info('survey realizations: %d : cell dmags' % i_realization)
#         z = np.random.normal(size=len(model.params.free))
#         zz = z / np.sqrt(D)
#         zz = fact.solve_Lt(zz)
#         zz = fact.apply_Pt(zz)
#         m = np.zeros(hp.nside2npix(nside))
#         m[:] = hp.UNSEEN
#         p = model.params.copy()
#         p.free = zz
#         m[dp.pixel_set] = p['mag'].full
#         ret.append(m)
#     return ret
