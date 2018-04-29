"""
Perform the observations, for a given cadence, in a given band
"""

from mx.DateTime import DateTimeFrom
import os

a = 2

# shared with ubercal studies
cadence_dir = '/data/nrl/LSST_GAIA_UBERCAL/cadence_db'

cadence_name, mjd_min, mjd_max, nside = get_input()
cadence_file = cadence_dir + os.sep + cadence_name + '.npy'
obs_file = get_data_fn('obs.npy')
nvisits_plot = get_data_fn('nvisits.png')

cmd = ['ubercal_observe.py', '-o', obs_file,
       '--mjd_min', str(mjd_min), '--mjd_max', str(mjd_max),
       '--nside', str(nside), 
       '--ncellsx', '1', '--ncellsy', '1', 
       '--fast', '1', '--norefpixs',
       cadence_file]
print ' '.join(cmd)
logged_subprocess(cmd)


seg_output = [(cadence_name, mjd_min, mjd_max, nside)]

