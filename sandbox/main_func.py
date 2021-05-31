from sub_func import *
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pickle import Pickler as pPk
from pickle import Unpickler as pUk

def Making_prep(folder, cadence):

    try:
        os.mkdir('Save_Metric')
    except:
        pass

    try:
        os.mkdir('Save_Metric/' + folder)
    except:
        pass

    try:
        os.mkdir('Save_Metric/' + folder + cadence)
    except:
        pass

    try:
        os.mkdir('Save_Metric/' + folder + cadence + '/fig')
    except:
        pass

    DATA = {'name': cadence}

    with open('Save_Metric/' + folder + cadence + '/save_data.dat', 'wb') as f:
        pPk(f).dump(DATA)

def Make_delay_metric(data, folder, cadence, bins=200, figsize=(16, 10)):

    """
    Give a statistique of cadence's gaps
    Give array for set of candence statistique

    Input :
    data (np.array) : data which containt cadence's data
    folder (str): folder which containt cadence
    cadence (str) : name of concerned cadence
    bins (int) : number of bins for histogramme in statistique
    figsize (tuple int) : tuple containt the figsize (a, b)*100 pixel

    """

    j = data['mjd']
    dj = j[1:] - j[:-1]
    dh = dj*24
    dm = dh*60

    plt.figure(figsize=(16, 10)).suptitle('Delay statistics')

    plt.subplot(231)
    plt.hist(dj, bins, color='r')
    plt.yscale('log')

    plt.subplot(234)
    plt.hist(dj[dj < 2], bins, color='r')
    plt.yscale('log')
    plt.xlabel('$\Delta$Mjd (day)')

    plt.subplot(232)
    plt.hist(dh[dh < 24], bins, color='r')
    plt.yscale('log')
    plt.title('Histogram $\Delta$Mjd of ' + cadence)

    plt.subplot(235)
    plt.hist(dh[dh < 5], bins, color='r')
    plt.yscale('log')
    plt.xlabel('$\Delta$Mjd (hour)')

    plt.subplot(233)
    plt.hist(dm[dh < 2], bins, color='r')
    plt.yscale('log')

    plt.subplot(236)
    plt.hist(dm[dm < 60], bins, color='r')
    plt.yscale('log')
    plt.xlabel('$\Delta$Mjd (min)')

    #plt.legend()
    plt.savefig('Save_Metric/' + folder + cadence + '/delay.png')

    data_delay = np.zeros(5).astype(float)

    data_delay[0] = np.sum(dj[dh < 1])
    data_delay[1] = np.sum(dj[dh < 9]) - np.sum(dj[dh < 1])
    data_delay[2] = np.sum(dj[dj < 1]) - np.sum(dj[dh < 9])
    data_delay[3] = np.sum(dj[dj < 7]) - np.sum(dj[dj < 1])
    data_delay[4] = np.sum(dj[dj >= 7])

    with open('Save_Metric/' + folder + cadence + '/save_data.dat', 'rb') as f:
        DATA = pUk(f).load()
    
    DATA['delay'] = data_delay

    with open('Save_Metric/' + folder + cadence + '/save_data.dat', 'wb') as f:
        pPk(f).dump(DATA)

def Make_pixel_metric(data, folder, cadence, argsMV={'cmap':'cool', 'cbar':False}, nside=64, recov=7, nb_day=30, fps=20, SUB=320, BAND=['gri', 'griz', 'g', 'r', 'i', 'z'], FMT=['.:k', '.-k', '.:r', '.:g', '.:b', '.:y'], figsize=(16, 10)):

    """
    Give metrics in rapport to activation of pixel in sky
        - Film of activ pixel under 7 days
        - %activ sky in function of days
        - give data for futur metrics (%activ sky and moy of activ day in row)

    data (np.array) : data cadence
    folder (str) : name fo the set of the cadence
    cadence (str) : name of the cadence
    argsMV (dict) : dict which containt args for Mollview
    nside (int) : size in pixel of the radius of sky in Mollview
    recov (int) : days until we put a pixel in activ (i dnt knw how to say ths)
    nb_day (int) : days when you want to applicate the metric
    fps (int) : number of wanted fps for he video
    SUB (int) : sub format for mollview
    BAND (list) : list of str with wanted band for metric
    FMT (list) : fmt wanted for metric
    figsize (tuple) : size of the figure (a, b)*100 pixel
    """

    argsMV['nest'] = True
    
    mjd_i = int(data['mjd'][0])
    mjd_f = int(data['mjd'][-1])

    t = np.arange(mjd_i, mjd_f)
    floor = np.floor(data['mjd']+5./24.)

    hpx0 = np.zeros(hp.nside2npix(nside)).astype(float) + recov
    tot_pxl = np.size(hpx0)
    f = FocPS()
    p = f.pixellize()

    HPX = []
    HPXs = []
    ADDpix = []
    t = t[:nb_day]

    MET = np.zeros((np.size(BAND), np.size(t)))

    for band in BAND:
        HPX.append(np.copy(hpx0))
        HPXs.append(np.copy(hpx0 * 0))
        ADDpix.append([])

    for k, ti in enumerate(t):
        I = np.where(floor == ti)
        d = data[I]

        plt.figure(figsize=figsize).suptitle('{} - [{}]'.format(cadence, ti))

        for i, band in enumerate(BAND):

            hpxN = make_hpmoll(d, HPX[i], recov, band, f, p, SUB + 1 + i, nside, argsMV)
            HPX[i] = np.copy(hpxN)

            MET[i, k] = np.size(np.where(HPX[i] != hp.UNSEEN)[0])

            fini = HPXs[i][hpxN == hp.UNSEEN]
            ADDpix[i] += list(fini[fini != 0])
            
            HPXs[i][hpxN == hp.UNSEEN] = 0
            HPXs[i][hpxN != hp.UNSEEN] += 1

        plt.savefig('Save_Metric/' + folder + cadence + '/fig/fig' + str(k) + '.png')
        plt.close()


    #Make the film with the figs
    path_folder = 'Save_Metric/' + folder + cadence + '/'
    create_film(nb_day, fps, path_folder+'fig/', prefixe='fig', extension='png')

    #Make duration activ pixel metric hsito
    plt.figure(figsize=figsize).suptitle('Count of duration of activ pixel')
    for i, band in enumerate(BAND):
        plt.subplot(SUB + 1 + i)
        plt.hist(ADDpix[i], 200, color='r')
        plt.title(band)

    plt.savefig('Save_Metric/' + folder + cadence + '/Metric_duration_activ_pixel.png')

    #Make the %sky activ pixel metric
    plt.figure(figsize=figsize)

    for fmt, band, met in zip(FMT, BAND, MET):

        plt.plot(t, met/tot_pxl*100, fmt, label=band)

    plt.xlabel('Day')
    plt.ylabel('% of activ pixel in sky')
    plt.title('Metric of Activ pixel of ' + cadence)
    plt.legend()

    plt.savefig(path_folder + 'pc_activ_sky.png')

    #Save data for cadence set metric
    moy_pcs = np.zeros(np.size(BAND))
    moy_act = np.zeros(np.size(BAND))

    for i, met in enumerate(MET):
        moy_pcs[i] = np.mean(met)/tot_pxl*100
        moy_act[i] = np.mean(ADDpix[i])

    with open(path_folder + 'save_data.dat', 'rb') as f:
        DATA = pUk(f).load()

    DATA['moy_pcs'] = moy_pcs
    DATA['moy_act'] = moy_act

    with open(path_folder + 'save_data.dat', 'wb') as f:
        pPk(f).dump(DATA)
    
def Make_repat_metric(data, folder, cadence):

    """

    """

    spy = 60*60*24 #sec per year
    slewTime = data['slewTime']

    #slew due at change filters
    dIt = np.zeros(np.size(data) - 1)

    for band in 'ugrizy':
        I = data['band'] == band
        dI = I[1:].astype(int) - I[:-1].astype(int)
        dI[dI == -1] = 0
        dIt += dI

    nbr_change_filter = np.sum(dIt)
    dIt = np.concatenate((np.array([0]), dIt))

    #Tout les temps :
    tota_t = find_total_time(data)
    deno = 100.0/spy/tota_t

    expt_t = np.sum(data['exptime'])*deno #Open Shutter time
    filt_t = np.sum(slewTime[dIt == 1])*deno #CHange filt time
    slew_t = np.sum(slewTime)*deno - filt_t #slew time
    read_t = np.sum(data['visitTime'])*deno - expt_t #readout time
    gaph_t = GiveGapHour(data)*deno #Gaps hour
    gapd_t = 100.0 - (expt_t + filt_t + slew_t + read_t + gaph_t) #Gaps day

    resume = [expt_t, slew_t, filt_t, read_t, gaph_t, gapd_t]

    with open('Save_Metric/' + folder + cadence + '/save_data.dat', 'rb') as f:
        DATA = pUk(f).load()

    DATA['pourcent'] = resume
    DATA['number_change_filter'] = nbr_change_filter
    DATA['number_mesure'] = np.size(data)

    with open('Save_Metric/' + folder + cadence + '/save_data.dat', 'wb') as f:
        pPk(f).dump(DATA)

    
#Metric pour les set de cadence

def make_set_metric(folder):

    listF = os.listdir('Save_Metric/' + folder)
    newL = []
    for f in listF:
        if ".png" in f:
            pass
        else:
            newL.append(f)
    listF = newL

    expt_name, expt_valr = np.zeros(np.size(listF)), np.zeros(np.size(listF))
    slew_name, slew_valr = np.zeros(np.size(listF)), np.zeros(np.size(listF))
    var_valr = np.zeros((np.size(listF), 6))

    delay = np.zeros((np.size(listF), 5))
    moyp_name, moyp_valr = np.zeros(np.size(listF)), np.zeros(np.size(listF))
    moyt_name, moyt_valr = np.zeros(np.size(listF)), np.zeros(np.size(listF))

    for i, fi in enumerate(listF):

        with open('Save_Metric/' + folder + fi + '/save_data.dat', 'rb') as f:
            DATA = pUk(f).load()

        var_valr[i] = DATA['pourcent']
        expt_valr[i] = DATA['pourcent'][0]
        slew_valr[i] = DATA['pourcent'][1]
        moyp_valr[i] = DATA['moy_pcs'][1]
        moyt_valr[i] = DATA['moy_act'][1]
        delay[i] = DATA['delay']

    name = np.array(listF)
    expt_name, slew_name, moyp_name, moyt_name = np.copy(name), np.copy(name), np.copy(name), np.copy(name)

    expt_index = np.argsort(expt_valr)
    slew_index = np.argsort(slew_valr)
    moyp_index = np.argsort(moyp_valr)
    moyt_index = np.argsort(moyt_valr)

    expt_name, expt_valr = expt_name[expt_index], expt_valr[expt_index]
    slew_name, slew_valr = slew_name[slew_index][::-1], slew_valr[slew_index][::-1]
    moyp_name, moyp_valr = moyp_name[moyp_index], moyp_valr[moyp_index]
    moyt_name, moyt_valr = moyt_name[moyt_index], moyt_valr[moyt_index]

    fig, ax = plt.subplots(figsize=(16, 10))
    plot = ax.plot(expt_valr, np.arange(np.size(expt_name)), c='r', ls=':', marker='.')
    ax.set_yticks(np.arange(np.size(expt_name)))
    ax.set_yticklabels(expt_name)
    ax.set_title('% Exptime')
    plt.savefig('Save_Metric/' + folder + 'Metric_Exptime.png')

    fig, ax = plt.subplots(figsize=(16, 10))
    plot = ax.plot(slew_valr, np.arange(np.size(slew_name)), c='r', ls=':', marker='.')
    ax.set_yticks(np.arange(np.size(slew_name)))
    ax.set_yticklabels(slew_name)
    ax.set_title('% Slewtime')
    plt.savefig('Save_Metric/' + folder + 'Metric_Slewtime.png')

    fig, ax = plt.subplots(figsize=(16, 10))
    plot = ax.plot(moyp_valr, np.arange(np.size(moyp_name)), c='r', ls=':', marker='.')
    ax.set_yticks(np.arange(np.size(moyp_name)))
    ax.set_yticklabels(moyp_name)
    ax.set_title('Mean %active pixel sky')
    plt.savefig('Save_Metric/' + folder + 'Metric_pc_active_sky.png')

    fig, ax = plt.subplots(figsize=(16, 10))
    plot = ax.plot(moyt_valr, np.arange(np.size(moyt_name)), c='r', ls=':', marker='.')
    ax.set_yticks(np.arange(np.size(moyt_name)))
    ax.set_yticklabels(moyt_name)
    ax.set_title('Mean duration of active pixel')
    plt.savefig('Save_Metric/' + folder + 'Metric_duration_pixel.png')
    
    #BAR
    categories = ['Open Shutter', 'Slew', 'Filter', 'Read', 'Gap in day', 'Gap > 1 day']
    labels = np.copy(expt_name)
    data = var_valr[expt_index]
    data_cum = data.cumsum(axis=1)
    cate_color = np.array(['r', 'orange', 'y', 'g', 'b', 'purple'])
    
    fig, ax = plt.subplots(figsize=(16, 10))
    ax.xaxis.set_visible(True)
    ax.set_xlim(0, np.sum(data, axis=1).max())
    ax.set_xlabel('% of Total Night time')

    for i, (colname, color) in enumerate(zip(categories, cate_color)):

        w = data[:, i]
        start = data_cum[:, i] -w
        rects = ax.barh(labels, w, left=start, height=0.5, label=colname, color=color)

        for i in range(1, 10):
            if i != 5:
                ax.axvline(i*10, color='k', ls=':', alpha=0.7)
            else:
                ax.axvline(i*10, color='k', ls='--', alpha=0.7)

        ax.legend(ncol=len(categories), bbox_to_anchor=(0, 1), loc='lower left', fontsize='small')
        plt.title('For ' + folder[:-1])

        plt.savefig('Save_Metric/' + folder + 'Metric_bartime.png')

    #DELAY
    cate = ['<1h', '<9h', '<1j', '<7j', '>7j']
    x = np.arange(len(cate))
    width = 0.06

    fig, ax = plt.subplots()
    RECT = []
    for i, d in enumerate(delay):

        rect = ax.bar(x + width/2 + i*width, d, width)
        RECT.append(rect)

    ax.set_title('Statistic of gaps')
    ax.set_xticks(x)
    ax.set_xticklabels(cate)
    #ax.legend()

    #for i, rect in enumerate(RECT):
    #    ax.bar_label(rect, padding=3)

    fig.tight_layout()
    plt.savefig('Save_Metric/' + folder + 'Metric_gap.png')

#Argparse

def parser_args_cadence():
    parser = argparse.ArgumentParser()
    parser.add_argument('data', type=np.ndarray)
    parser.add_argument('folder', type=str)
    parser.add_argument('cadence', type=str)
    return parser.parse_args()

def parser_args_set():
    parser = argparse.ArgumentParser()
    parser.add_argument('folder', type=str)
    parser.add_argument('cadence', type=str)
    return parser.parse_args()
    
    
#Main cadence

def main_cadence(data, folder, cadence):
    Making_prep(folder, cadence)
    Make_delay_metric(data, folder, cadence)
    Make_pixel_metric(data, folder, cadence, nb_day=16)
    Make_repat_metric(data, folder, cadence)

#Main set of cadence
    
def main_set(folder, cadence):
    make_set_metric(folder, cadence)
