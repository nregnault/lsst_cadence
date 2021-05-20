"""
for each cadence, collapse all frames into a movie
"""

import os
import os.path as op

cadence_name, mjd_min, mjd_max, nside, sn = get_input()

obs_frame_0 = glob_parent('00000.png')[0]
frame_dir = op.dirname(obs_frame_0)
movie_fn = get_data_fn('%s_%d_%d_%d_%s.mp4' % (cadence_name, mjd_min, mjd_max, nside, sn))

cmd = ['ffmpeg-hi10-heaac', '-y', '-i', frame_dir + os.sep + '%05d.png',
       '-pix_fmt', 'yuv420p',
       '-c:v', 'libx264', '-movflags', '+faststart',
       movie_fn]
logged_subprocess(cmd)

seg_output = [(cadence_name, mjd_min, mjd_max, nside, sn)]

