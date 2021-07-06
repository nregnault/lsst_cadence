#multiplex cross_prod group_by 'global_metric[-1]'

import os
import numpy as np
import matplotlib.pyplot as plt
from pickle import Unpickler as pUk
from ALF_cadfunc.class_html import MakeHtml as MH

change = 'trivial 12'

plt.switch_backend('Agg')
plt.ioff()

num = get_input()

sets_dir = glob_parent('*')
sup_dir = get_data_fn('.')

#print('\n###')
#for seti in sets_dir:
#    print(seti, end='\n\n')

try:
    os.mkdir(sup_dir + '/save_metric')
except:
    pass

COLOR_vr = {'fbs_1_3':'#FF0000',
           'fbs_1_4':'#CC0044',
           'fbs_1_5':'#990099',
           'fbs_1_6':'#4400CC',
            'fbs_1_6_1':'#220066',
           'fbs_1_7':'#0000FF',
            'fbs_1_7_1':'#000099',
           'altsched_new':'#00FF00',
           'altsched_orig':'#009900'}


#MEAN GRI SUPER GLOBAL
plt.figure(figsize=(16, 10))

for set_file in sets_dir:

    with open(set_file + '/save_metric/save_data.dat', 'rb') as f:
        supDATA = pUk(f).load()

    version = supDATA['version']
    duration = supDATA['duration_mean']
    meansky = supDATA['meansky']
    alpha_facecolor = supDATA['alpha_facecolor']

    color_i = COLOR_vr[version]
    facecolor = []
    for alpha in alpha_facecolor:
        facecolor.append(color_i + alpha)

    plt.scatter(duration, meansky, edgecolor=color_i, facecolor=facecolor)

for ver in COLOR_vr.keys():
    plt.plot([] ,[], ls='', marker='.', color=COLOR_vr[ver], label=ver)

plt.scatter([], [], edgecolor='k', facecolor='#00000000', label='Rolling')

plt.legend()
plt.xlabel('Mean activating duration (days)')
plt.ylabel('Activation area (%sky)')
plt.savefig(sup_dir + '/save_metric/Super_GRI_mean.png')

#MEDIAN GRI SUPER GLOBAL
plt.figure(figsize=(16, 10))

for set_file in sets_dir:

    with open(set_file + '/save_metric/save_data.dat', 'rb') as f:
        supDATA = pUk(f).load()

    version = supDATA['version']
    duration = supDATA['duration_medi']
    meansky = supDATA['meansky']

    plt.scatter(duration, meansky, color=COLOR_vr[version])

for ver in COLOR_vr.keys():
    plt.plot([] ,[], ls='', marker='.', color=COLOR_vr[ver], label=ver)

plt.legend()
plt.xlabel('Mean activating duration (days)')
plt.ylabel('Activation area (%sky)')
plt.savefig(sup_dir + '/save_metric/Super_GRI_median.png')

#make html

text_hover = ['Home']
link_hover = [sup_dir + '/save_metric/super_global.html']

text_nav_global = []
link_nav_global = []

for fold in sets_dir:

    with open(fold + '/save_metric/save_data.dat', 'rb') as f:
        supDATA = pUk(f).load()

    version = supDATA['version']

    text_hover = ['Global']
    link_hover = [sup_dir + '/save_metric/super_global.html']

    text_nav_global.append(version)
    link_nav_global.append(fold + '/save_metric/' + version + '.html')

    text_nav = []
    link_nav = []

    for cad in os.listdir(fold + '/save_metric/cadence_html/'):

        text_nav.append(cad[:-5])
        link_nav.append(fold + '/save_metric/cadence_html/' + cad)

    html = MH(fold + '/save_metric/' + version)
    html.add_Hover(text_hover, link_hover)
    html.add_Navigation(text_nav, link_nav, 'Cadences :')

    html.begin_body()
    html.add_label('Metrics of {} :'.format(version))

    html.add_label('Bartime metric : ', 'h2')
    html.add_pict(fold + '/save_metric/Metric_bartime.png')

    html.add_label('Slewtime metric : ', 'h2')
    html.add_pict(fold + '/save_metric/Metric_Slewtime.png')

    html.add_label('Exptime metric : ', 'h2')
    html.add_pict(fold + '/save_metric/Metric_Exptime.png')

    html.add_label('Graphic with mean activation pixel duration (for GRI) : ', 'h2')
    html.add_pict(fold + '/save_metric/Metric_gri_Meanduration.png')

    html.add_label('Graphic with median activation pixel duration (for GRI) :', 'h2')
    html.add_pict(fold + '/save_metric/Metric_gri_Medianduration.png')
    
    html.add_label('Gap statistique : ', 'h2')
    html.add_pict(fold + '/save_metric/Metric_gap.png')

    html.close()



html = MH(sup_dir + '/save_metric/super_global')
html.add_Hover(text_hover, link_hover)
html.add_Navigation(text_nav_global, link_nav_global, 'Version : ')
html.begin_body()
html.add_label('Global Metric')

html.add_label('Mean duration', 'h2')
html.add_pict(sup_dir + '/save_metric/Super_GRI_mean.png')

html.add_label('Median duration', 'h2')
html.add_pict(sup_dir + '/save_metric/Super_GRI_mean.png')

html.close()
