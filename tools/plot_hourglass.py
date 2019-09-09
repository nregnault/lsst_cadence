#!/usr/bin/env python


import os
import os.path as op
import sys

import numpy as np
import pylab as pl


colors = {
    'u': 'b',
    'g': 'g',
    'r': 'r',
    'i': 'orange',
    'z': 'darkred',
    'y': 'k'
    }


colors = {
    'u': 'navy',
    'g': 'cyan',
    'r': 'b',
    'i': 'y',
    'z': 'r',
    'y': 'k'
    }


def hourglass_plot(d, mjd_lims=(61400, 61700), savefig=None):
    fig = pl.figure()
    for b in "ugrizy":
        idx = d['band'] == b
        dd = d[idx]
        mjd = np.floor(dd['mjd']+5./24.)
        pl.plot(mjd, dd['mjd']-mjd, color=colors[b],
                marker='.', ls='', markersize=0.5)
        pl.plot(mjd[0], dd['mjd'][0]-mjd[0], color=colors[b],
                marker='.', ls='', label=b)
    pl.xlabel('MJD [days]')
    pl.legend(loc='upper center', ncol=6,
              bbox_to_anchor=(0.5, 1.1))
    ax = pl.gca()
    #    ax.get_yaxis().set_visible(False)
    ax.get_yaxis().set_ticks([])
    pl.ylabel('exposures')
    if mjd_lims is not None:
        pl.xlim(mjd_lims)
    if savefig is not None:
        fig.savefig(savefig, bbox_inches='tight')


def main(filename, savefig=None, mjd_lims=None):
    cad = np.load(filename)
    hourglass_plot(cad, savefig=savefig, mjd_lims=mjd_lims)
