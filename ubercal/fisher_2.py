#multiplex cross_prod where 'observe[1] in ubercalfits[0]' group_by '(observe[0], observe[2], observe[3], tuple(ubercalfits[0]), ubercalfits[1])'

"""
Build fisher matrix (for the moment, just for the default model)
and factorize it. 
"""

import os
import os.path as op



cadence_name, band, nside, yrs, fluxrefs = get_input()


obs_file = glob_parent('obs.npy.npz')
plot_dir = op.dirname(get_data_fn('toto'))
nb_realizations = 100

print obs_file

cmd = ['ubercal_fisher.py', 
       '--output_dir', plot_dir,
       '--dump_cls', 
       '--dump_maps',
       #       '--dump_inv_fisher', plot_dir,
       '--realizations', str(nb_realizations), 
       '--plot_dir', plot_dir, 
       '--refpixs', str(fluxrefs),
       '--fit', '1'] + obs_file
print ' '.join(cmd)
logged_subprocess(cmd)

seg_output = [(cadence_name, band, nside, yrs, fluxrefs)]

