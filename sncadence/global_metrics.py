#multiplex cross_prod group_by 'observe[:-1]'

"""
global metrics for the cadence
"""

import os
import os.path as op

a = 3

cadence_name, mjd_min, mjd_max, nside = get_input()

obs_file = glob_parent('obs.npy.npz')[0]
plot_dir = op.dirname(get_data_fn('toto'))

cmd = ['global_cadence_stats.py', 
       '-O', plot_dir, 
       '--nside', str(nside), 
       obs_file]
print ' '.join(cmd)
logged_subprocess(cmd)

seg_output = [(cadence_name, mjd_min, mjd_max, nside)]
