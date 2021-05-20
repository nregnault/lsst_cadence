#!/usr/bin/env python

"""
Download a subset of the GAIA DR2
borrowed from Marc Betoule/ 
"""

from astropy.utils.data import download_file, clear_download_cache, download_files_in_parallel
import os
import os.path as op
import shutil
import re
import logging
logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',
                    level=logging.INFO)
from argparse import ArgumentParser

import subprocess
import numpy as np
from croaks import convert
import healpy as hp
# from imageproc.reference_catalog import SkyCat




parse_url = re.compile('href="(.*)">')
gaia_archive = 'http://cdn.gea.esac.esa.int/Gaia/gdr2/gaia_source/csv/'
directory = '../data/'

#    from saunerie.astrotools import sexagecimal_to_decimal as s2d
# fields = {'P330E': (247.86, 30.26),
#           'P177D': (239.78, 47.74),
#           'GD71': (88.00, 16.00),
#           'SNAP 2': (244.87, 55.68),
#           'GD153': (194.26, 22.00),
#           'Grid-1': (s2d('06:30:00.00',hourly=True),s2d('+14:20:00.0')),
#           'Grid-2': (s2d('21:00:00.00',hourly=True),s2d('+10:00:00.0')),
#           'D1'    : (s2d('02:26:00.00',hourly=True),s2d('-04:30:00.0')),
#           'D2'    : (s2d('10:00:28.60',hourly=True),s2d('+02:12:21.0')),
#           'D3'    : (s2d('14:19:28.01',hourly=True),s2d('+52:40:41.0')),
#           'D4'    : (s2d('22:15:31.67',hourly=True),s2d('-17:44:05.7')),
# }
    

def get_pixlist(level=7, max_declination=30., plot=False):
    """
    """
    nside = hp.order2nside(level)
    npix = hp.nside2npix(nside)
    ipix = np.arange(npix)
    theta,phi = hp.pix2ang(nside, ipix, nest=True)
    #    dec = -theta * 180. / np.pi + 90.
    dec = -np.degrees(theta - np.pi / 2.)
    idx = dec < 35.

    if plot:
        m = np.zeros(npix)
        m[idx] = 1.
        hp.mollview(m, nest=True)
    return ipix[idx]
    
def index_file_to_url_list(index_file):
    with open(index_file) as f:
        l = f.readlines()
    return [gaia_archive+parse_url.findall(_l)[0] for _l in l if ('href' in _l and '.gz' in _l)]

def down(url, dest):
    print 'downloading %s' % url
    fname = download_file(url, cache=False)
    shutil.move(fname, dest)

def readcsvgz(fname):
    content = subprocess.check_output(['zcat', fname]).splitlines()
    names = content[0].split(',')
    vals = [[convert(v) for v in _l.split(',')] for _l in content[1:]]  
    return np.rec.fromrecords(vals, names=names)

def readsub(fname):
    from io import StringIO
    def cfloat(s):
        return float(s.strip() or 'nan')

    def pmerror(s):
        return float(s.strip() or 0) * 1e3

    def mag(s):
        return float(s.strip() or 0) *1e3

    def magerror(s):
        return 1.0857362047581294 / float(s.strip() or np.inf) * 1e3

    def parallax(s):
        return float(s.strip() or 0) * 10
    
    output_format = [('ra', 'f8', cfloat, None), # deg
                     ('dec', 'f8', cfloat, None), # deg
                     ('ra_error', 'u2', pmerror, None), # mas -> micro as
                     ('dec_error', 'u2', pmerror, None), # mas -> micro as
                     ('pmra', 'f4', cfloat, None), # mas / an
                     ('pmdec', 'f4', cfloat, None), # mas / an
                     ('parallax', 'i2', parallax, None), # mas -> 0.1mas
                     ('pmra_error', 'u2', pmerror, None), # mas/an -> micro as / an
                     ('pmdec_error', 'u2', pmerror, None), # mas/an -> micro as / an
                     ('phot_g_mean_mag', 'u2', mag, 'g'), # mag -> mmag
                     ('phot_bp_mean_mag', 'u2', mag, 'bp'), # mag -> mmag
                     ('phot_rp_mean_mag', 'u2', mag, 'rp'), # mag -> mmag
                     ('phot_g_mean_flux_over_error', 'u2', magerror, 'g_error'), # mmag
                     ('phot_bp_mean_flux_over_error', 'u2', magerror, 'bp_error'), # mmag
                     ('phot_rp_mean_flux_over_error', 'u2', magerror, 'rp_error') # mmag
    ]
    
    content = StringIO(unicode(subprocess.check_output(['zcat', fname])))
    l = content.readline()
    names = l.split(',')
    cols = [names.index(c[0]) for c in output_format]
    return np.loadtxt(content, dtype={'names': [e[3] if e[3] else e[0] for e in output_format],
                                'formats': [e[1] for e in output_format]},
                      delimiter=',', usecols=cols,
                      converters = dict(zip(cols, [e[2] for e in output_format])))

def url_to_healpix(urllist, level=12):
    ids = np.array([[int(v) for v in h.split('.')[-3].split('_')[2:]] for h in urllist])    
    div = 2**35 * 4**(12 - level)
    return np.rec.fromarrays([ids[:, 0] / div, ids[:, 1] / div], names=['start','stop'])

def interesting_files(urllist, pixlist, level=12):
    healpix = url_to_healpix(urllist, level=level)
    return np.array([np.any([start <= h and h <= stop for h in pixlist]) for start, stop in healpix])
    
def get_pixlist_for_fields(fields, radius=2, level=12):
    import healpy
    pixlist = []
    radius = np.radians(radius)
    nside = healpy.order2nside(level)
    for ra, dec in fields:
        vec = healpy.ang2vec(ra, dec, lonlat=True)
        pixlist.append(healpy.query_disc(nside, vec, radius, inclusive=True, fact=4, nest=True))
    pixlist = np.hstack(pixlist)
    return np.unique(pixlist)


def main(output_dir='./'):
    # retrieve index 
    logging.info('downloading index from: %s' % gaia_archive)
    index_file = download_file(gaia_archive, cache=True)
    urllist = index_file_to_url_list(index_file)

    # file selection 
    logging.info('%d files' % len(urllist))    
    healpix = url_to_healpix(urllist)
    level=7    
    pixlist = get_pixlist(level=4)
    selection = interesting_files(urllist, pixlist, level=4)
    nfiles = selection.sum()
    logging.info('%d interesting files' % nfiles)

    if not op.isdir(output_dir):
        os.makedirs(output_dir)
        
    # download
    for _u,url in enumerate(np.array(urllist)[selection], 1):
        print 'downloading file %d/%d' % (_u, nfiles)
        down(url, 'temp.gz')
        r = readsub('temp.gz')
        print r.__class__
        np.save(output_dir + os.sep + 'c_%d.npy' % _u, r)


if __name__ == "__main__":
    parser = ArgumentParser(description='download a subset of GAIA DR2')
    parser.add_argument('-O', '--output_dir', default='./',
                        help='output directory')
    args = parser.parse_args()
    
    main(output_dir=args.output_dir)
    
    
    #        skycat.add_safe(readsub('temp.gz'))
    #    skycat  = SkyCat(outdir='../data/gaiadr2_healpix')
    #toto = readsub('./GaiaSource_1000172165251650944_1000424567594791808.csv.gz')
    #toto.tofile('toto.npy')
    #urls = index_file_to_url_list(index_file)
    #for u in urls:
    #    down(u, u.replace(vizir_ftp, directory))
