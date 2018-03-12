"""
Perform the observations, for a given cadence, in a given band
"""

cadence_dir = ''

cadence_name, mjd_min, mjd_max, band, nside = get_input()
# log = get_log(cadence_name, mjd_min, mjd_max, band, cadence_dir=cadence_dir)
obs_file = get_data_fn('obs.npy')

cmd = ['ubercal_observe', '-o', obs_file,
       '--start', mjd_min, '--stop', mjd_max,
       '--band', band, cadence_dir]


output_seg = [(cadence_name, mjd_min, mjd_max, band, nside]

