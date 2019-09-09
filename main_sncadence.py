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

from mx.DateTime import DateTimeFrom

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
            'altsched_good_weather', 
            'altsched_rolling_good_weather', 
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

# Additional project cadences
cadences += """mothra_2049 kraken_2042 kraken_2044 nexus_2097""".split()

# Altsched + twilight 
cadences += """alt_sched_twi""".split()

# P. Yoachim experiments
# cadences += """blobs_same_10yrs""".split()
cadences += """blobs_mix_zmask10yrs blobs_same_10yrs blobs_same_10yrs blobs_same_zmask10yrs rolling_10yrs 
cadence_mix_10yrs rolling_mix_10yrs rolling_mix_75_10yrs
tight_mask_10yrs tight_mask_simple_10yrs tms_drive_10yrs tms_roll_10yrs
""".split()

# P. Gris new altscheds 
cadences += """altsched_new_ddf1-90.0_18.0_30.0_20_1_10yrs altsched_new_ddf1-90.0_18.0_30.0_20_2_10yrs altsched_new_ddf1-90.0_3.0_30.0_20_1_10yrs altsched_new_ddf1-90.0_3.0_30.0_20_1_10yrs_rolling altsched_new_ddf1-90.0_3.0_30.0_20_2_10yrs altsched_new_ddf1-90.0_3.0_30.0_20_2_10yrs_rolling altsched_new_ddf1-90.0_3.0_30.0_30_2_10yrs""".split()

# cadences = """altsched_good_weather kraken_2026""".split()

cadences += "new_wfd1yrs".split()

cadences = ['baseline_1exp_nopairs_10yrs', 'baseline_1exp_pairsame_10yrs', 'baseline_1exp_pairsmix_10yrs', 'baseline_2exp_pairsame_10yrs', 'baseline_2exp_pairsmix_10yrs', 
            'ddf_0.23deg_1exp_pairsmix_10yrs', 'ddf_0.70deg_1exp_pairsmix_10yrs', 'ddf_pn_0.23deg_1exp_pairsmix_10yrs', 'ddf_pn_0.70deg_1exp_pairsmix_10yrs', 
            'exptime_1exp_pairsmix_10yrs', 
            'baseline10yrs', 
            'big_sky10yrs', 'big_sky_nouiy10yrs', 
            'gp_heavy10yrs', 
            'newA10yrs', 'newB10yrs', 
            'roll_mod2_mixed_10yrs', 'roll_mod3_mixed_10yrs', 'roll_mod6_mixed_10yrs', 
            'simple_roll_mod10_mixed_10yrs', 'simple_roll_mod2_mixed_10yrs', 'simple_roll_mod3_mixed_10yrs', 'simple_roll_mod5_mixed_10yrs', 
            'twilight_1s10yrs', 'altsched_1exp_pairsmix_10yrs', 'rotator_1exp_pairsmix_10yrs', 
            'hyak_baseline_1exp_nopairs_10yrs', 'hyak_baseline_1exp_pairsame_10yrs']

cadences = ['pontus_2573', 'rolling_10yrs_opsim', 'rolling_mix_10yrs_opsim']

cadences = ['very_alt2_rm5illum20_10yrs', 'very_alt2_rm5illum40_10yrs', 'very_alt3_rm5illum20_10yrs',
            'very_alt3_rm5illum40_10yrs', 'very_alt10yrs', 'very_alt2_rm5illum25_10yrs',
            'very_alt2_rm5illum50_10yrs', 'very_alt3_rm5illum25_10yrs', 'very_alt3_rm5illum50_10yrs',
            'very_alt2_rm5illum15_10yrs', 'very_alt2_rm5illum30_10yrs', 'very_alt3_rm5illum15_10yrs',
            'very_alt3_rm5illum30_10yrs', 'very_alt_rm510yrs', 'noddf_1exp_pairsame_10yrs',
            'desc_ddf_pn_0.70deg_1exp_pairsmix_10yrs', 'fc1exp_pairsmix_ilim30_10yrs', 'fc1exp_pairsmix_ilim60_10yrs',
            'fc1exp_pairsmix_ilim15_10yrs', 'stuck_rolling10yrs', 'shortt_2ns_1ext_pairsmix_10yrs',
            'shortt_2ns_5ext_pairsmix_10yrs', 'shortt_5ns_5ext_pairsmix_10yrs', 'shortt_5ns_1ext_pairsmix_10yrs',
            'simple_roll_mod2_mixed_10yrs', 'roll_mod2_sdf0.2mixed_10yrs', 'simple_roll_mod3_sdf0.2mixed_10yrs',
            'roll_mod2_sdf0.1mixed_10yrs', 'roll_mod3_sdf0.2mixed_10yrs', 'roll_mod3_sdf0.1mixed_10yrs',
            'simple_roll_mod5_sdf0.2mixed_10yrs', 'roll_mod6_sdf0.2mixed_10yrs', 'roll_mod6_sdf0.1mixed_10yrs',
            'simple_roll_mod10_sdf0.2mixed_10yrs', 'roll_mod2_sdf0.10mixed_10yrs', 'roll_mod2_sdf0.05mixed_10yrs',
            'simple_roll_mod2_sdf0.20mixed_10yrs', 'roll_mod3_sdf0.05mixed_10yrs', 'roll_mod2_sdf0.20mixed_10yrs',
            'roll_mod3_sdf0.20mixed_10yrs', 'simple_roll_mod3_sdf0.20mixed_10yrs', 'roll_mod3_sdf0.10mixed_10yrs',
            'roll_mod6_sdf0.05mixed_10yrs', 'roll_mod6_sdf0.20mixed_10yrs', 'roll_mod6_sdf0.10mixed_10yrs',
            'simple_roll_mod10_sdf0.20mixed_10yrs']

# all cadences, including the latests ones from June/July 2019
# with the more convervative weather files
cadences = ['altsched_18_-90_30', 'big_sky10yrs',
            'kraken_2035', 'roll_mod6_sdf0.05mixed_10yrs',
            'tms_roll_10yrs', 'altsched_18_-90_40',
            'big_sky_nouiy0yrs', 'kraken_2036',
            'roll_mod6_sdf0.10mixed_10yrs', 'too_pairsmix_rate100_10yrs',
            'altsched_1exp_pairsmix_10yrs', 'big_sky_nouiy10yrs',
            'kraken_2042', 'roll_mod6_sdf0.20mixed_10yrs',
            'too_pairsmix_rate10_10yrs', 'altsched_good_weather',
            'blobs_mix_zmask10yrs', 'kraken_2044',
            'rotator_1exp_pairsmix_10yrs', 'too_pairsmix_rate1_10yrs',
            'altsched_new_ddf1-90.0_18.0_30.0_20_1_10yrs',
            'blobs_same_10yrs', 'minion_1016',
            'shortt_2ns_1ext_pairsmix_10yrs', 'too_pairsmix_rate50_10yrs',
            'altsched_new_ddf1-90.0_18.0_30.0_20_2_10yrs',
            'blobs_same_zmask10yrs', 'mothra_2045',
            'shortt_2ns_5ext_pairsmix_10yrs', 'twilight_1s10yrs',
            'altsched_new_ddf1-90.0_3.0_30.0_20_1_10yrs',
            'bluer_footprint10yrs', 'mothra_2049',
            'shortt_5ns_1ext_pairsmix_10yrs', 'very_alt10yrs',
            'altsched_new_ddf1-90.0_3.0_30.0_20_1_10yrs_rolling',
            'cadence_mix_10yrs', 'newA10yrs',
            'shortt_5ns_5ext_pairsmix_10yrs',
            'very_alt2_rm5illum15_10yrs',
            'altsched_new_ddf1-90.0_3.0_30.0_20_2_10yrs', 'colossus_2664',
            'newB10yrs', 'simple_roll_mod10_sdf0.20mixed_10yrs',
            'very_alt2_rm5illum20_10yrs',
            'altsched_new_ddf1-90.0_3.0_30.0_20_2_10yrs_rolling',
            'colossus_2665', 'new_wfd1yrs',
            'simple_roll_mod2_sdf0.20mixed_10yrs',
            'very_alt2_rm5illum25_10yrs',
            'altsched_new_ddf1-90.0_3.0_30.0_30_2_10yrs', 'colossus_2667',
            'nexus_2097', 'simple_roll_mod3_sdf0.20mixed_10yrs',
            'very_alt2_rm5illum30_10yrs', 'alt_sched',
            'ddf_0.23deg_1exp_pairsmix_10yrs',
            'noddf_1exp_pairsame_10yrs',
            'simple_roll_mod5_sdf0.20mixed_10yrs',
            'very_alt2_rm5illum40_10yrs', 'altsched_rolling_good_weather',
            'ddf_0.70deg_1exp_pairsmix_10yrs', 'pontus_2002',
            'stability_-10offset_42seed10yrs',
            'very_alt2_rm5illum50_10yrs', 'alt_sched_rolling',
            'ddf_pn_0.23deg_1exp_pairsmix_10yrs', 'pontus_2489',
            'stability_10offset_42seed10yrs',
            'very_alt3_rm5illum15_10yrs', 'alt_sched_twi',
            'ddf_pn_0.70deg_1exp_pairsmix_10yrs', 'pontus_2502',
            'stability_180offset_42seed10yrs',
            'very_alt3_rm5illum20_10yrs', 'baseline0yrs',
            'dec_1exp_pairsmix_10yrs', 'pontus_2573',
            'stability_1offset_42seed10yrs', 'very_alt3_rm5illum25_10yrs',
            'baseline10yrs', 'desc_ddf_pn_0.70deg_1exp_pairsmix_10yrs',
            'presto_10yrs', 'stability_1offset_43seed10yrs',
            'very_alt3_rm5illum30_10yrs', 
            'exptime_1exp_pairsmix_10yrs', 'presto_third_10yrs',
            'stability_1offset_44seed10yrs', 'very_alt3_rm5illum40_10yrs',
            'baseline_1exp_nopairs_10yrs', 'fc1exp_pairsmix_ilim15_10yrs',
            'rolling_10yrs', 'stability_30offset_42seed10yrs',
            'very_alt3_rm5illum50_10yrs', 
            'fc1exp_pairsmix_ilim30_10yrs', 'rolling_10yrs_opsim',
            'stability_365offset_42seed10yrs', 'very_alt_rm510yrs',
            'baseline_1exp_pairsame_10yrs',
            'fc1exp_pairsmix_ilim60_10yrs', 'rolling_mix_10yrs',
            'stuck_rolling10yrs', 'weather_0.10c_10yrs',
            'feature_baseline_10yrs',
            'rolling_mix_10yrs_opsim',
            'templates_w_1.0_1exp_pairsmix_10yrs', 'weather_0.20c_10yrs',
            'baseline_1exp_pairsmix_10yrs',
            'feature_rolling_half_mask_10yrs', 'rolling_mix_75_10yrs',
            'templates_w_2.0_1exp_pairsmix_10yrs', 'weather_0.30c_10yrs',
            'baseline2018a', 'feature_rolling_twoThird_10yrs',
            'roll_mod2_sdf0.05mixed_10yrs',
            'templates_w_3.0_1exp_pairsmix_10yrs', 'weather_0.40c_10yrs',
            'gp_heavy0yrs',
            'roll_mod2_sdf0.10mixed_10yrs',
            'templates_w_4.0_1exp_pairsmix_10yrs', 'weather_0.60c_10yrs',
            'baseline_2exp_pairsame_10yrs', 'gp_heavy10yrs',
            'roll_mod2_sdf0.20mixed_10yrs',
            'templates_w_5.0_1exp_pairsmix_10yrs', 'weather_0.70c_10yrs',
            'hyak_baseline_1exp_nopairs_10yrs',
            'roll_mod3_sdf0.05mixed_10yrs', 'tight_mask_10yrs',
            'weather_0.80c_10yrs', 'baseline_2exp_pairsmix_10yrs',
            'hyak_baseline_1exp_pairsame_10yrs',
            'roll_mod3_sdf0.10mixed_10yrs', 'tight_mask_simple_10yrs',
            'weather_0.90c_10yrs', 'big_sky0yrs', 'kraken_2026',
            'roll_mod3_sdf0.20mixed_10yrs', 'tms_drive_10yrs',
            'weather_1.10c_10yrs',]



def get_tasks(cadences, nside):
    ret = []
    mjd_min = DateTimeFrom('2022-01-01').mjd
    mjd_max = DateTimeFrom('2032-12-31').mjd
    seasons = [(mjd_min, mjd_max)]
    
    for c in cadences:
        for begin, end in seasons:
            ret.extend([(c, begin, end, ns) for ns in nside])
    return ret


def start_server(P, address, debug=1):
    filename = 'pipe_sncadence.pkl'
    import cPickle
    with open(filename, 'w') as f:
        cPickle.dump(P, f)
    
    cmd = ['pipeletd', '-n', '-l', str(debug), 
           '-a', address[0], '-p', str(address[1]), 
           filename]
    print ' '.join(cmd)
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
    tasks = get_tasks(cadences, nside=args.nside)
    P.push(observe=tasks)
    
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
        W,t = launch_process(P, args.process, 
                             log_level=log_level, 
                             address=('', 56001))
        W.run()
        
    return P

if __name__ == '__main__':
    P = main()
    

