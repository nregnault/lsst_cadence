"""
analyze the cadence
"""

import os
import os.path as op

cadence_name, mjd_min, mjd_max, nside, sn = get_input()

obs_file = glob_parent('obs.npy.npz')[0]
plot_dir = op.dirname(get_data_fn('toto'))
ebv_mask = '/sps/lsst/cadence/LSST_SN_CADENCE/static/ebv_mask_%d.npy' % nside
tabulated_nsn = '/sps/lsst/cadence/LSST_SN_CADENCE/static/tabulated_nsn_%d.npz' % nside
lc_templates = {'faint': '/sps/lsst/cadence/LSST_SN_CADENCE/static/lc_templates.npy',
                'normal': '/sps/lsst/cadence/LSST_SN_CADENCE/static/lc_templates_normal.npy'}
vmax_nsn = {'faint':  '15',
            'normal': '25'}


cmd = ['animate_cadence.py', '-O', plot_dir, 
       '--ebv-mask', ebv_mask, obs_file, 
       '--nsn', tabulated_nsn,
       '--min_cadence', '0.25',
       '--vmax_nsn', vmax_nsn[sn], 
       '--nside', str(nside),
       '--lc_template', lc_templates[sn]]
print ' '.join(cmd)
logged_subprocess(cmd)

seg_output = [(cadence_name, mjd_min, mjd_max, nside, sn)]
