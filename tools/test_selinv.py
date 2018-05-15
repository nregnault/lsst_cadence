"""
"""

import os
import numpy as np
import pylab as pl

import logging
logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

import scipy.sparse
import scipy.sparse.linalg
from sksparse import cholmod
from saunerie import selinv


def generate_posdef(n=1000, nv=5000):
    # generate a matrix that is (supposedly) posdef
    i = np.random.randint(n, size=nv)
    j = np.random.randint(n, size=nv)
    v = np.random.uniform(0., 1., size=nv)
    H = scipy.sparse.coo_matrix((v,(i,j)), shape=(n,n))
    D = scipy.sparse.dia_matrix((1.E-6 * np.ones(n), [0]), shape=(n,n))
    M = (D + H.T * H).tocsc()
    return M


def generate_realization(fact, n):
    D = fact.D()
    z = np.random.normal(size=n)
    zz = z / np.sqrt(D)
    zz = fact.solve_Lt(zz)
    zz = fact.apply_Pt(zz)
    return zz


def factorize_and_invert(M):
    n,_ = M.shape
    
    # factorize it
    logging.info('cholesky')
    fact = cholmod.cholesky(M)
    
    # invert it
    logging.info('brute force invert')    
    inv_dense = np.linalg.inv(M.todense())
    logging.info('sparse invert')        
    inv_sparse = scipy.sparse.linalg.inv(M).todense()
    
    # inverse of the diagonal
    logging.info('selinv')
    r = selinv.selinv(M, P=fact.P())

    # MC inverse
    logging.info('MC inverse')
    rlz = np.vstack([generate_realization(fact, n) for i in xrange(50000)])
    
    return fact, r, np.diagonal(inv_dense), np.diagonal(inv_sparse), rlz


def main():
    M = generate_posdef(n=1000, nv=5000)
    fact, r, inv_dense, inv_sparse, rlz = factorize_and_invert(M)
    v = np.var(rlz, axis=0)
    return M, fact, r, inv_dense, inv_sparse, v, rlz
    
