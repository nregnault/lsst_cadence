#multiplex cross_prod group_by 'cadence_metric[-1]'

from ALF_cadfunc.main_function import global_metric as gm
from ALF_cadfunc.class_html import MakeHtml as MH
import os
from pickle import Unpickler as pUk

change = 'trivial changement 10'

version = get_input()

if '.' in version:
    version = '_'.join(version.split('.'))

cads_dir = glob_parent('')
set_dir = get_data_fn('.')

print('\n###')
for cad in cads_dir:
    print(cad, end='\n\n')

try:
    os.mkdir(set_dir + '/save_metric')
except:
    pass

gm(set_dir, cads_dir, version)

#make html

text_nav = []
link_nav = []

for fold in cads_dir:

    with open(fold + '/save_metric/save_data.dat', 'rb') as f:
        DATA = pUk(f).load()

    text_nav.append(DATA['name'])
    link_nav.append(DATA['name'] + '.html')

try:
    os.mkdir(set_dir + '/save_metric/cadence_html')
except:
    pass

for i, fold in enumerate(cads_dir):

    directory = fold + '/save_metric/'
    
    text_hover = [version]
    link_hover = [set_dir + '/save_metric/' + version + '.html']

    html = MH(set_dir + '/save_metric/cadence_html/' + text_nav[i])
    html.add_Hover(text_hover, link_hover)
    html.add_Navigation(text_nav, link_nav, 'Cadences :')
    html.begin_body()
    html.add_label('Metrics of {} :'.format(text_nav[i]))

    html.add_label('Film of activ pixel :', 'h2')
    html.add_mp4(directory + 'film.mp4')

    html.add_label('Pourcentage of active sky in GRI :', 'h2')
    html.add_pict(directory + 'pc_activ_sky_gri.png')

    html.add_label('Statistic of duration pixel :', 'h2')
    html.add_pict(directory + 'Metric_duration_activ_pixel.png')

    html.add_label('Statistic of gaps :', 'h2')
    html.add_pict(directory + 'delay.png')

    html.close()


set_output = [(version, 0)]
