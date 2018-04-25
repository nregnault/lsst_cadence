#!/usr/bin/env python

import os
import os.path as op
import sys

import numpy as np
import pylab as pl
import argparse
import logging
logging.basicConfig(fmt="%(asctime)s %(levelname)s %(message)",
                    level=logging.INFO)
from saunerie import salt2, snsim
from pycosmo.cosmo import CosmoLambda

cosmo = CosmoLambda()


def get_pulse_shapes(bands,
                     X1=-2., Color=0.2, DayMax=0.,                     
                     restframe_phase_range=(-20., 45.),
                     redshifts=np.arange(0.05, 1.2, 0.01),
                     filename='salt2.npz'):
    """
    """
    lcmodel = snsim.init_lc_model(bands, filename=filename)
    pmin, pmax = restframe_phase_range

    l = []
    for z in redshifts:
        mjd_min = np.floor(pmin * (1+z) + DayMax)
        mjd_max = np.ceil(pmax * (1+z) + DayMax)
        mjd = np.arange(mjd_min, mjd_max+1, 1.)
        sn = np.rec.fromrecords([(z, cosmo.dL(z), X1, Color, DayMax)],
                                names=['z', 'dL', 'X1', 'Color', 'DayMax'])
        n = len(mjd)
        o = np.ones(n)
        for bn in bands:
            b = [bn] * n
            val = lcmodel(sn, mjd, b)
            l.append(np.rec.fromarrays((mjd, val, z*o, b),
                                       names=['mjd', 'val', 'z', 'band']))
            
    return np.hstack(l)
    


if __name__ == "__main__":
    bands = ['LSSTPG::' + b for b in "grizy"]
    r = get_pulse_shapes(bands)
    np.save('lc_templates.npy', r)
    
