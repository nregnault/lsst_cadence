"""
Perform the observations, for a given cadence, in a given band
"""

from mx.DateTime import DateTimeFrom
import os

a = 2
#comment='fbs v1.5 and v1.6 cadences (for Alex & Satya)'
comment='fbs v1.7'
expose(['comment'])


# shared with ubercal studies
cadence_dir = '/sps/lsst/cadence/LSST_SN_CADENCE/cadence_db'

cadence_name, mjd_min, mjd_max, nside = get_input()
cadence_file = cadence_dir + os.sep + cadence_name + '.npy'
obs_file = get_data_fn('obs.npy')
nvisits_plot = get_data_fn('nvisits.png')

cmd = ['ubercal_observe.py', '-o', obs_file,
       '--mjd_min', str(mjd_min), '--mjd_max', str(mjd_max),
       '--nside', str(nside), 
       '--ncellsx', '1', '--ncellsy', '1', 
       '--fast', '1',
       cadence_file]
print ' '.join(cmd)
logged_subprocess(cmd)

seg_output = [(cadence_name, mjd_min, mjd_max, nside, sn) for sn in ['faint', 'normal']]

