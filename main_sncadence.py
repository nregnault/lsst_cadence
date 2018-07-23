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
observe -> analyze;
observe -> global_metrics;
"""

log_level = logging.INFO
code_dir = op.abspath('./sncadence')
prefix = '/sps/lsst/cadence/LSST_SN_CADENCE/'
sql_file = prefix + os.sep + '.sqlstatus'

# Jan 2018 cadences
cadences = ['alt_sched',   
            'alt_sched_rolling',
            'feature_baseline_10yrs',
            'feature_rolling_half_mask_10yrs',
            'feature_rolling_twoThird_10yrs',
            'minion_1016']

# white paper cadences
cadences += """baseline2018a colossus_2667 kraken_2035
pontus_2489 colossus_2664 mothra_2045 colossus_2665
kraken_2026 pontus_2002 kraken_2036 pontus_2502""".split()

# AltSched variants
cadences += """altsched_18_-90_30 altsched_18_-90_40""".split()

def get_tasks(cadences, nside):
    ret = []
    mjd_min = DateTimeFrom('2022-01-01').mjd
    mjd_max = DateTimeFrom('2032-12-31').mjd
    seasons = [(mjd_min, mjd_max)]
    
    for c in cadences:
        for begin, end in seasons:
            ret.append((c, begin, end, nside))
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
                        dest='nside', type=int, default=64)
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
    tasks = get_tasks(cadences, nside=args.nside)
    P.push(observe=tasks)
    
    if args.debug:
        W,t = launch_interactive(P, log_level=log_level)
        W.run()
    else:
        W,t = launch_process(P, args.process, 
                             log_level=log_level, 
                             address=('', 56001))
        W.run()
        
    return P

if __name__ == '__main__':
    P = main()
    

