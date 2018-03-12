#!/usr/bin/env python 

"""
"""

import pipelet.pipeline as pipeline
from pipelet.launchers import launch_interactive, launch_process
import os
import os.path as op
import sys
import shutil

from mx.DateTime import DateTimeFrom

import matplotlib
matplotlib.use('Agg')
import pylab as pl 
pl.interactive(0)


import logging
import numpy as np
from croaks import NTuple


pipedot="""
observe;
"""

log_level = logging.INFO
code_dir = op.abspath('./ubercal')
prefix = '/data/nrl/LSST_GAIA_UBERCAL/'
sql_file = prefix + os.sep + '.sqlstatus'

cadences = ['alt_sched',   
            'alt_sched_rolling',
            'feature_baseline_10yrs',
            'feature_rolling_half_mask_10yrs',
            'feature_rolling_twoThird_10yrs',
            'minion_1016']
mjds = [DateTimeFrom(str(yr) + '-01-01').mjd for yr in xrange(2022, 2033)]


def get_tasks(cadences, mjds, bands, nside):
    ret = []
    seasons = np.repeat(mjds, 2)[1:-1].reshape((-1,2))
    print seasons
    
    for c in cadences:
        for begin, end in seasons:
            for b in bands:
                ret.append((c, begin, end, b, nside))
    return ret


def main():
    """
    run the pipeline
    """
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--debug', 
                        help='start jobs in interactive mode',
                        dest='debug', action='store_true', default=False)
    parser.add_argument('--nside', 
                        help='healpy nside',
                        dest='nside', type=int, default=1024)
    parser.add_argument('-p', '--process', metavar='N', 
                        help='launch jobs as local parallel processes',
                        type=int, default=1)
    parser.add_argument('-D', '--dump_pipeline',
                        help='dump the pipeline structure in a dot file',
                        action='store_true', dest='dump', default=False)
    args = parser.parse_args()


    # pipeline instance 
    P = pipeline.Pipeline(pipedot, 
                          code_dir=code_dir,
                          prefix=prefix, 
                          sqlfile=sql_file)
    
    # print out the pipeline structure
    if args.dump:
        P.to_dot('pipeline.dot')
        sys.exit(0)
        
    #    tasks = [('minion_1016', 59580., 59945., 'r', 1024)]
    tasks = get_tasks(cadences, mjds, ['i', 'z'], nside=args.nside)
    P.push(observe=tasks)
    
    if args.debug:
        W,t = launch_interactive(P, log_level=log_level)
        W.run()
    else:
        W,t = launch_process(P, args.process, 
                             log_level=log_level, 
                             address=('', 56000))
        W.run()
        
    return P

if __name__ == '__main__':
    P = main()
    

