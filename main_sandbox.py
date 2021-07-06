#!/usr/bin/env python 

"""
"""

import pipelet.pipeline as pipeline
from pipelet.launchers import launch_interactive, launch_process, launch_pbs, launch_ccage_worker_2
import os
import os.path as op
import sys
import shutil
import subprocess

#from mx.DateTime import DateTimeFrom

import matplotlib
matplotlib.use('Agg')
import pylab as pl 
pl.interactive(0)


import logging
import numpy as np
from croaks import NTuple

PROJECT='lsst'
CPU_TIME = "160:00:00"
MEMORY_SIZE = "16G"
SCRATCH_SIZE = "30G"
JOB_DIR = "/sps/lsst/cadence/LSST_SN_CADENCE/bqs/scripts"
LOG_DIR = "/sps/lsst/cadence/LSST_SN_CADENCE/bqs/logs"
PLATFORM_SPECS = ""
QUEUE = 'mc_longlasting'

#dot l
pipedot="""
cadence_metric -> global_metric -> super_global_metric;
"""

log_level = logging.INFO
code_dir = op.abspath('./Test0')
#prefix = '/data/lsst/cadence/Test0/'
prefix = '/data/lsst/cadence/prods/result_pipeline/'
sql_file = prefix + os.sep + '.sqlstatus'

racine_version = '/data/lsst/cadence/cadence_db/'

version = ['fbs_1.3',
           'fbs_1.4',
           'fbs_1.5',
           'fbs_1.6',
           'fbs_1.6.1',
           'fbs_1.7',
           'fbs_1.7.1',
           'altsched_orig',
           'altsched_new']

#version = ['fbs_1.7.1']

cadences = []

for ver in version:
    listcad = os.listdir(racine_version + ver)
    for f in listcad:
        if '.npy' in f:
            cadences.append(ver + '/' + f)
            #print(cadences[-1])

"""

cadences = [
    'altsched_new/altsched_new_ddf1-90.0_18.0_30.0_20_1_10yrs.npy',
    'altsched_new/altsched_new_ddf1-90.0_3.0_30.0_20_1_10yrs_rolling.npy',
    'altsched_new/altsched_new_ddf1-90.0_3.0_30.0_30_2_10yrs.npy',
    'fbs_1.7/rolling_scale0.2_nslice2_v1.7_10yrs.npy',
    'fbs_1.7/ddf_dither0.00_v1.7_10yrs.npy',
    'fbs_1.7/footprint_0_v1.710yrs.npy']

"""


#def get_tasks(cadences, nside):
#    ret = []
#    mjd_min = DateTimeFrom('2022-01-01').mjd # was 2022-01-01
#    mjd_max = DateTimeFrom('2032-12-31').mjd # was 2032-12-31
#    seasons = [(mjd_min, mjd_max)]
#    for c in cadences:
#        for begin, end in seasons:
#            ret.extend([(c, begin, end, ns) for ns in nside])
#    return ret


#ignore
def start_server(P, address, debug=1):
    filename = 'pipe_sncadence.pkl'
    import cPickle
    with open(filename, 'w') as f:
        cPickle.dump(P, f)
    
    cmd = ['pipeletd', '-n', '-l', str(debug), 
           '-a', address[0], '-p', str(address[1]), 
           filename]
    print(' '.join(cmd))
    #    subprocess.Popen(cmd).communicate()[0]

    

def main():
    """
    run the pipeline
    """
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--debug', 
                        help='start jobs in interactive mode',
                        dest='debug', action='store_true', default=False)
    parser.add_argument('--nside', nargs='+',
                        help='healpy nside',
                        dest='nside', type=int, default=[64,])
    parser.add_argument('-p', '--process', metavar='N', 
                        help='launch jobs as local parallel processes',
                        type=int, default=1)
    parser.add_argument('-D', '--dump_pipeline',
                        help='dump the pipeline structure in a dot file',
                        action='store_true', dest='dump', default=False)
    parser.add_argument('-w', '--add-workers', metavar='N',
                        help='submit N additional jobs without launching a new server', 
                        type=int)
    parser.add_argument('-S', '--start-server', default=False, 
                        action='store_true',
                        help='start the pipelet server (typically on ccwsge)')
    
    parser.add_argument('--project', default="lsst", type=str,
                        help='project to which these jobs belong')
    parser.add_argument('--port', default=50010,
                        type=int,
                        help='port the scheduler is listening on')
    parser.add_argument('--cpu', default="72:00:00", type=str,
                        help='CPU time soft limit (batch)')
    parser.add_argument('--vmem', default="15G", type=str,
                        help='virtual memory limit (batch)')
    parser.add_argument('--scratch', default="25G", type=str,
                        help='scratch size soft limit (batch)')
    parser.add_argument('--mc', default=None, type=int,
                        help='multicore option (mandatory if submitting to mc queue)')
    parser.add_argument('--queue', default=None, # longlasting
                        help='queue name')
    parser.add_argument('-N', '--ntasks_per_worker', default=None, type=int,
                        help='maximum number of tasks per worker')
    
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
    #    tasks = get_tasks(cadences, nside=args.nside)
    tasks = [(cn, -1, cn.split('/')[0]) for cn in cadences]
    P.push(cadence_metric=tasks)

    #tasks = [(st) for st in sets]
    #P.push(main_set=tasks)
    
    if args.debug:
        W,t = launch_interactive(P, log_level=log_level)
        W.run()
    elif args.start_server:
        start_server(P, address=('ccwsge1348.in2p3.fr', args.port))
    elif args.add_workers:
        launch_ccage_worker_2(P, args.add_workers, 
                              address=('ccwsge1348.in2p3.fr', args.port), 
                              project=args.project,
                              job_name='lsst_cadence',
                              # job_dir=JOB_DIR,
                              # log_dir=LOG_DIR,
                              queue=args.queue,
                              cpu_time=args.cpu, 
                              vmem=args.vmem, 
                              scratch=args.scratch,
                              multicores=args.mc,
                              log_level=log_level,
                              ntasks_per_worker=args.ntasks_per_worker)
    else:
        print(' @@@ ', P, args.process, ('', args.port))
        W,t = launch_process(P, args.process, 
                             log_level=log_level, 
                             address=('', args.port))
        W.run()
        
    return P

if __name__ == '__main__':
    P = main()
    

