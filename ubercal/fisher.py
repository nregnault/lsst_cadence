"""
Build fisher matrix (for the moment, just for the default model)
and factorize it. 
"""

import os
import os.path as op

a = 6

cadence_name, mjd_min, mjd_max, band, nside = get_input()

obs_file = glob_parent('obs.npy.npz')[0]
plot_dir = op.dirname(get_data_fn('toto'))
nb_realizations = 100

cmd = ['ubercal_fisher.py', 
       '--output_dir', plot_dir,
       '--dump_cls', 
       '--dump_maps',
       #       '--dump_inv_fisher', plot_dir,
       '--realizations', str(nb_realizations), 
       '--plot_dir', plot_dir, 
       obs_file]
print ' '.join(cmd)
logged_subprocess(cmd)

seg_output = [(cadence_name, mjd_min, mjd_max, band, nside)]

