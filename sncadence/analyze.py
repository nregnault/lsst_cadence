"""
analyze the cadence
"""

import os
import os.path as op

a = 5

cadence_name, mjd_min, mjd_max, nside = get_input()

obs_file = glob_parent('obs.npy.npz')[0]
plot_dir = op.dirname(get_data_fn('toto'))
ebv_mask = '/data/nrl/LSST_SN_CADENCE/static/ebv_mask.npy'
tabulated_nsn = '/data/nrl/LSST_SN_CADENCE/static/tabulated_nsn.npz'
lc_templates = '/data/nrl/LSST_SN_CADENCE/static/lc_templates.npy'

cmd = ['animate_cadence.py', '-O', plot_dir, 
       '--ebv-mask', ebv_mask, obs_file, 
       '--nsn', tabulated_nsn,
       '--min_cadence', '0.25',
       '--vmax_nsn', '30', 
       '--lc_template', lc_templates]
print ' '.join(cmd)
logged_subprocess(cmd)

seg_output = [(cadence_name, mjd_min, mjd_max, nside)]
