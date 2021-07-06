import os
import numpy as np
from .lsst import FocalPlaneSimple as FocPS
import healpy as hp
from PIL import Image
from scipy.interpolate import interp1d as ipd

def change_name(cadence0):
    list_sup = ['.npy', '_10yrs', '10yrs']

    cadence = cadence0

    for st in list_sup:
        cadence = ''.join(cadence.split(st))

    return cadence

def make_hpmoll(d, hpx0, recov, band, f, p, sub, nside, argsMV):

    hpx = np.copy(hpx0)
    hpx[hpx0 != hp.UNSEEN] = hpx0[hpx0 != hp.UNSEEN] + 1
    hpx[hpx > recov] = hp.UNSEEN

    ra, dec = np.array([]), np.array([])

    for b in band:
        ra  = np.concatenate((ra,  d[d['band'] == b]['Ra']))
        dec = np.concatenate((dec, d[d['band'] == b]['Dec']))

    if len(ra) > 0:

        for r, d in zip(ra, dec):
            cells = p.copy()
            f.to_uv(cells, (r, d))

            for c in cells:
                poly = np.vstack([c['p0'], c['p1'], c['p2'], c['p3']])

                try:
                    ipix = hp.query_polygon(nside, poly, nest=True, inclusive=False)
                    hpx[ipix] = 0
                except:
                    continue

    argsMV['title'] = band
    argsMV['sub'] = sub

    hp.mollview(hpx, min=0, max=recov, **argsMV)

    return hpx

def create_film(nb_frame, fps, folder, prefixe, extension, init0):

    spf = 1.0/fps #sec per frame
    mspf = int(spf*1000) #msec per frame
    path_in = folder + prefixe
    path_out = folder  + 'AnimPIL.gif'

    framename = []

    for i in range(nb_frame):

        framename.append(path_in + str(init0 + i) + '.' + extension)

    img, *imgs = [Image.open(f) for f in framename] #opening des images avec PIL
    img.save(fp=path_out, format='GIF', append_images=imgs, save_all=True, duration=mspf, loop=0) #save gif !

def find_total_time(data):

    """
    Fit a hourglass plot and give the total observable time
    """

    mjd = data['mjd']

    floor = np.floor(mjd + 5./24.)
    expo = mjd - floor

    delta_up = np.append(floor[1:] - floor[:-1], 1)
    derde_up = np.where(delta_up == 1)
    x_up = expo[derde_up]

    delta_dw = floor[1:] - floor[:-1]
    derde_dw = np.where(delta_dw == 1)[0] + 1
    derde_dw = np.concatenate((np.array([0]), derde_dw))
    x_dw = expo[derde_dw]

    floor_x = floor[derde_dw]

    #extract of wrong point
    I_up = np.where(x_up > 0.33)
    I_dw = np.where(x_dw < 0.06)

    floor_up, x_up = floor_x[I_up], x_up[I_up]
    floor_dw, x_dw = floor_x[I_dw], x_dw[I_dw]

    #Interpollation
    func_up = ipd(floor_up, x_up)
    func_dw = ipd(floor_dw, x_dw)

    t = np.arange(floor_x[0], floor_x[-1])
    delta = func_up(t) - func_dw(t)

    return np.sum(delta)

def GiveGapHour(data):

    mjd = data['mjd']
    dj = mjd[1:] - mjd[:-1]
    dh = dj*24
    dm = dh*60
    ds = dm*60

    Ii = np.where(dm > 5)
    Is = np.where(dh[Ii] < 10)

    return np.sum(ds[Ii][Is])
