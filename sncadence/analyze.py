"""
analyze the cadence
"""

import os
import os.path as op

cadence_name, mjd_min, mjd_max, nside, sn = get_input()

zmax = 0.6
#comment='new reprocessing with new altsched cadences'
#expose(['zmax', 'comment'])
#comment='reprocessing for Alex & Satya (2019-2020 cadence paper)'
comment='processing of the fbs 1.7 cadences [zmax=0.6]'
expose(['zmax', 'comment'])

# sn = 'normal'

obs_file = glob_parent('obs.npy.npz')[0]
plot_dir = op.dirname(get_data_fn('toto'))
ebv_mask = '/sps/lsst/cadence/LSST_SN_CADENCE/static/ebv_mask_%d.npy' % nside
tabulated_nsn = '/sps/lsst/cadence/LSST_SN_CADENCE/static/tabulated_nsn_%d.npz' % nside
lc_templates = {'faint': '/sps/lsst/cadence/LSST_SN_CADENCE/static/lc_templates.npy',
                'normal': '/sps/lsst/cadence/LSST_SN_CADENCE/static/lc_templates_normal.npy'}
vmax_nsn = {'faint':  '10',  # 15
            'normal': '15'}  # 25


cmd = ['animate_cadence.py', '-O', plot_dir, 
       '--ebv-mask', ebv_mask, obs_file, 
       '--nsn', tabulated_nsn,
       '--min_cadence', '0.25', # was 0.01
       '--vmax_nsn', vmax_nsn[sn], 
       '--nside', str(nside),
       '--zmax', str(zmax),
       '--lc_template', lc_templates[sn],
       '--cadence-name', cadence_name,
       #       '--drop-snr-req',
       #       '--drop-early-late-req',]
       #       '--drop-min-cadence-req',
   ]
print ' '.join(cmd)
logged_subprocess(cmd)

seg_output = [(cadence_name, mjd_min, mjd_max, nside, sn)]
