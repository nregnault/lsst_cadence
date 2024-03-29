"""
Perform the observations, for a given cadence, in a given band
"""

import os


cadence_dir = '/sps/lsst/cadence/UBERCAL/cadence_db' 

cadence_name, yr, band, nside = get_input()
cadence_file = cadence_dir + os.sep + cadence_name + '.npy'
obs_file = get_data_fn('obs.npy')
nvisits_plot = get_data_fn('nvisits.png')

cmd = ['ubercal_observe.py', '-o', obs_file,
       #       '--mjd_min', str(mjd_min), '--mjd_max', str(mjd_max),
       '--band', band, '--nside', str(nside), 
       '--ncellsx', '1', '--ncellsy', '1', 
       '--plots', nvisits_plot,
       '--year', str(yr),
       cadence_file]
    
print ' '.join(cmd)
logged_subprocess(cmd)


seg_output = [(cadence_name, yr, band, nside)]

