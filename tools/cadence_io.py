#!/usr/bin/env python 

"""
This is the fragile part of the whole thing, when we parse the output
of the pipeline and build a consistent a simple directory tree
containing the results.

This part of the code will probably change often. 
"""


import os
import os.path as op
import sys

import re
import glob

import numpy as np

import logging
logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',
                    level=logging.INFO)


from argparse import ArgumentParser 

from ubercal.cadence import Cadence, SummaryData




class FileNotFound(Exception): pass


def pipeletize_path(p):
    return p.replace('.', '_')


def strip_suffix(fn):
    i = fn.rfind('.')
    if i>=0:
        return fn[:i]


    
def _load_summary_file(root_dir, sntype, cadences=None):
    """
    Load the summary file and return a dictionary of cadences 
    (not fully filled yet, since the useful data (plots, movie ...)
    are going to be loaded by the _load_data function (below)
    
    .. note :: all the information is already in the movie/global_metrics segment. 
               not sure the summary files are that important anymore. 
    """
    
    # can be called several times with the same cadence dictionary
    if cadences is None:
        cadences = {}

    # at the moment, the pipeline spits out something called final_state
    fn = glob.glob(root_dir + os.sep + '*_' + sntype + os.sep + 'final_state_sim_history.npy')
    if len(fn) == 0:
        raise FileNotFound('unable to find the summary file')
    elif len(fn) > 1:
        raise FileNotFound('more than one file matches the search pattern: %s' % pattern)
    data = np.load(fn[0])    

    # match the cadence names specified in the summary ntuple
    # in order to extract the short cadence name and the FBS version
    # a side effect of this is that is removes all the <10 years simulations.
    matches = [re.match('(.+)_?(v\d*\.\d*)_10yrs', cn.decode()) for cn in data['cadence_name']]
    for i,r in enumerate(matches):
        if r is None:
            logging.warning('unable to parse cadence name: %s' % data[i]['cadence_name'])
            continue
        name, version = r.groups()
        name = name.strip('_ ')
    
        c = cadences.get(name, None)
        if c is None:
            c = Cadence(name, version)
            cadences[name] = c
        c.summary_data[sntype] = SummaryData(i,data)
        
    return cadences


def _load_data(root_dir, cadences):
    """
    """
    def _loadit(root_dir, cad, suffix, sntype=None):
        """
        """
        pattern = pipeletize_path('*%s*%s*' % (c.name, c.version_tag))
        if sntype is not None:
            pattern += sntype
        pattern = root_dir + os.sep + pattern + os.sep + '*' + suffix
        g = glob.glob(pattern)
        if len(g) > 0:
            if sntype is not None:
                g = list(filter(lambda x: ((sntype in x) and (re.search(c.name + '_?' + c.version_tag, x) is not None)), g))
            else:
                g = list(filter(lambda x: re.search(c.name + '_?' + c.version_tag, x) is not None, g))
            return g
            #                k = strip_suffix(op.basename(fn))
            #                v = op.relpath(fn)
            #                if k in target:
            #                    logging.warning('[_loadit] warning: key %s already in target' % k)
            #                target[k] = v
        else:
            raise FileNotFound('unable to locate global metric plots for cadence: %s, %s, %s' % (c.name, c.version_tag, pattern))
        
        return None
    
    lst = cadences if isinstance(cadences, list) else cadences.values()
    for c in lst:
        g = _loadit(root_dir, c, '.png', sntype='normal')
        k = [strip_suffix(op.basename(fn)) for fn in g]
        #        v = [op.relpath(fn) for fn in g]
        v = [op.abspath(fn) for fn in g]
        c.plots = dict(zip(k,v))
        
        g = _loadit(root_dir, c, '.npz', sntype='normal')
        k = [strip_suffix(op.basename(fn)) for fn in g]
        #        v = [op.relpath(fn) for fn in g]
        v = [op.abspath(fn) for fn in g]
        c.data = dict(zip(k,v))
        
        g = _loadit(root_dir, c, suffix='mp4')
        c.movie_file = {}
        for fn in g:
            if 'faint' in fn:
                c.movie_file['faint'] = op.abspath(fn)
            elif 'normal' in fn:
                c.movie_file['normal'] = op.abspath(fn)
            else:
                raise FileNotFound('do not know what to do with: %s' % fn)
            

def load_cadences(summary_data_dir, merged_data_dir):
    """
    load the cadences along with 
    """
    cadences = {}
    
    _load_summary_file(summary_data_dir, 'faint', cadences)
    _load_summary_file(summary_data_dir, 'normal', cadences)

    _load_data(merged_data_dir, cadences)

    return cadences


def dump_cadences(cads, root_dir, full_names=False):
    """
    """
    if not op.isdir(root_dir):
        os.makedirs(root_dir)

    import copy
    ret = copy.deepcopy(cads)
        
    for k,v in ret.items():
        logging.info('copying cadence: %s' %v.name)
        if full_names:
            cadence_dir = root_dir + os.sep + v.full_name
        else:
            cadence_dir = root_dir + os.sep + v.name

        for sntype in ['normal', 'faint']:
            target_dir = cadence_dir + os.sep + sntype
            if not op.isdir(target_dir):
                logging.info('creating: %s' % target_dir)
                os.makedirs(target_dir)

            tgt = target_dir + os.sep + op.basename(v.movie_file[sntype])
            if not op.isfile(tgt):
                os.link(op.realpath(v.movie_file[sntype]), tgt)
            v.movie_file[sntype] = tgt
            for k,p in v.plots.items():
                tgt = target_dir + os.sep + op.basename(p)
                if not op.isfile(tgt):
                    os.link(op.realpath(p), tgt)
                v.plots[k] = tgt
                
            for k,p in v.data.items():
                tgt = target_dir + os.sep + op.basename(p)
                if not op.isfile(tgt):
                    os.link(op.realpath(p), tgt)
                v.data[k] = tgt
            for sn,s in v.summary_data.items():
                tgt = target_dir + os.sep + 'summary_' + sn + '.npy'
                if not op.isfile(tgt):
                    np.save(target_dir + os.sep + 'summary_' + sn + '.npy', s.data)

    return ret



if __name__ == '__main__':
    parser = ArgumentParser(description='cadence_io.py: buid a pickled summary database for the cadences')
    parser.add_argument('--summary_dir', 
                        help='summary data directory')
    parser.add_argument('--merged_results_dir', 
                        help='directory that containts all the plots etc')
    parser.add_argument('-o', '--output', 
                        default='cadences.pkl',
                        help='pickled cadence db')
    parser.add_argument('-O', '--output-dir', 
                        default='./',
                        help='where to dump pickled cadence_db and plots')

    args = parser.parse_args()
    
    #    cads = load_cadences(summary_data_dir='summary_XKDZ7SI/data/', merged_data_dir='movies_3UZZB6Q/merge_results_HEQVO2Y/data/')
    cads = load_cadences(summary_data_dir=args.summary_dir, merged_data_dir=args.merged_results_dir)
    ret = dump_cadences(cads, args.output_dir + os.sep + 'plots')
    
    import pickle
    with open(args.output_dir + os.sep + 'cads_debug.pkl', 'wb') as f:
        pickle.dump(cads, f)
    with open(args.output_dir + os.sep + args.output, 'wb') as f:
        pickle.dump(ret, f)
        
    
