from ALF_cadfunc.main_function import cadence_metric as cm
import os
import numpy as np

#changement trivial pour reboot
change = 'changement trivial 15'

cadence_dir = '/data/lsst/cadence/cadence_db'
cadence_name, nb_day, version = get_input()
cadence_file = cadence_dir + os.sep + cadence_name

cadence_dir_s = get_data_fn('.')

cad = cadence_name.split('/')
folder = cad[0] + '/'
cadence = cad[1]

cm(np.load(cadence_file), cadence_dir_s, cadence, nb_day=nb_day)


frame_dir = cadence_dir_s + '/save_metric/fig/'
savem_dir = cadence_dir_s + '/save_metric/'
#movie_fn = get_data_fn('film.mp4')
picture = os.listdir(frame_dir)
picture.sort()
first_num = str(picture[0][3:8])

cmd = ['ffmpeg', '-y', '-start_number', first_num, '-i',
       frame_dir + 'fig%05d.png',
       '-pix_fmt', 'yuv420p',
       '-c:v', 'libx264', '-movflags', '+faststart',
       savem_dir + 'film.mp4']

logged_subprocess(cmd)

seg_output = [(cadence_name, version)]





