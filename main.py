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
observe -> fisher;
"""
# observe -> fisher;

log_level = logging.INFO
code_dir = op.abspath('./ubercal')
prefix = '/sps/lsst/cadence/UBERCAL' 
sql_file = prefix + os.sep + '.sqlstatus'

cadences = ['alt_sched',   
            'alt_sched_rolling',
            'feature_baseline_10yrs',
            'feature_rolling_half_mask_10yrs',
            'feature_rolling_twoThird_10yrs',
            'minion_1016']

# white paer cadences
cadences += """baseline2018a colossus_2667 kraken_2035
pontus_2489 colossus_2664 mothra_2045 colossus_2665
kraken_2026 pontus_2002 kraken_2036 pontus_2502""".split()

# AltSched variants
cadences += """altsched_18_-90_30 altsched_18_-90_40""".split()

# Additional project cadences
cadences += """mothra_2049 kraken_2042 kraken_2044 nexus_2097""".split()

cadences += """blobs_mix_zmask10yrs blobs_same_10yrs blobs_same_10yrs blobs_same_zmask10yrs rolling_10yrs 
cadence_mix_10yrs rolling_mix_10yrs rolling_mix_75_10yrs
tight_mask_10yrs tight_mask_simple_10yrs tms_drive_10yrs tms_roll_10yrs
""".split()


def get_tasks_mjd(cadences, mjd, bands, nside, standards):
    mjds = [DateTimeFrom(str(yr) + '-01-01').mjd for yr in xrange(2022, 2033, 1)]
    seasons = np.repeat(mjds, 2)[1:-1].reshape((-1,2))
    print seasons
    
    ret = []    
    for c in cadences:
        for begin, end in seasons:
            for b in bands:
                ret.append((c, begin, end, b, nside, standards))
    return ret


def get_tasks(cadences, bands, nside, nb_years=10):
    ret = []
    for c in cadences:
        for yr in np.arange(nb_years):
            for b in bands:
                ret.append((c,yr,b,nside))
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

    # get the task list 
    tasks = get_tasks(cadences, ['z'], nside=args.nside)
    P.push(observe=tasks)
    
    if args.debug:
        W,t = launch_interactive(P, log_level=log_level)
        W.run()
    else:
        W,t = launch_process(P, args.process, 
                             log_level=log_level, 
                             address=('', 56003))
        W.run()
        
    return P

if __name__ == '__main__':
    P = main()
    

