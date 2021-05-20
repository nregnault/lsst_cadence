#!/usr/bin/env python 


"""
"""

import os
import os.path as op
import sys
import re

import logging
logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',
                    level=logging.INFO)

import numpy as np
import pylab as pl

from adjustText import adjust_text

from argparse import ArgumentParser



def plot_nsn_vs_z(fn, vmin=0.3, vmax=0.5, cmap=pl.cm.jet, title='', savefig=None, plot_labels=True, **kw):
    """
    """
    d = np.load(fn)
    pl.figure(figsize=(16,9))
    #    db = load_db(fn=cadence_db)

    pl.plot(d['zmax'], d['nsn_tot'], ls='', color='gray', marker='.', alpha=0.5)
    texts, labels = [], []

    # plot the cadences once
    pl.plot(d['zmax'], d['nsn_tot'], 'k,')

    # scatter plot for
    cadence_names, rolling = [], []
    for c in d:
        cadence_names.append(c['cadence_name'].decode())
        rolling = ['rolling' in nm for nm in cadence_names]
    rolling = np.array(rolling).astype(bool)
    
    pl.scatter(d[~rolling]['zmax'], d[~rolling]['nsn_tot'], c=d[~rolling]['cadence'], vmin=vmin, vmax=vmax, cmap=cmap, label='non rolling')
    pl.scatter(d[rolling]['zmax'], d[rolling]['nsn_tot'], c=d[rolling]['cadence'], vmin=vmin, vmax=vmax, cmap=cmap, marker='^', label='rolling')
    pl.colorbar()
    pl.xlabel('$z_{max}$')
    pl.ylabel('$N_{SNe}$ [passing minimal sampling requirements]')
    pl.title(title)

    if 'xmin' in kw or 'xmax' in kw:
        pl.xlim((kw.get('xmin', None), kw.get('xmax', None)))
    if 'ymin' in kw or 'ymax' in kw:        
        pl.ylim((kw.get('ymin', None), kw.get('ymax', None)))

    
    if plot_labels == True:
        for c in d:
            #            if (c['zmax']>0.28 and c['nsn_tot']>120000.) or (c['nsn_tot']<100000):
            name = c['cadence_name'].decode()
            r = re.match('(.+)(_?v1.7).*10yrs', name)
            if r is None:
                logging.warning('unable to parse: %s' % name)
                continue
            else:
                nm, ver = r.groups()
            texts.append(pl.text(c['zmax'], c['nsn_tot'], nm, fontsize=12, va='center', ha='center', color='k'))
                
        adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red'))
    pl.grid(1)

    if savefig is not None:
        pl.gcf().savefig(savefig, bbox_inches='tight')

        
if __name__ == '__main__':
    """
    """
    parser = ArgumentParser(description='final summary plots')
    parser.add_argument('-o', '--output', default='nsn_z.png',
                        help='output file')
    parser.add_argument('--nolabels', 
                        default=False, action='store_true', 
                        help='do not label the strategies on the plot')
    parser.add_argument('--title', 
                        default='$N_{SN} vs z$ [normal SN]', 
                        help='specify figure title')
    parser.add_argument('histfile', 
                        help='input history file from which summary plots are made')
    
    args = parser.parse_args()
    pl.switch_backend('Agg')

    print 'LABELS: ', 
    plot_nsn_vs_z(args.histfile, plot_labels=not args.nolabels, title=args.title, savefig=args.output)

    #    plot_nsn_vs_z('summary_XKDZ7SI/data/6596_normal/final_state_sim_history.npy', plot_labels=False, title='$N_{SN} vs z$ [normal SN]', savefig='nsn_z_normal_nolabels.png')
    #    plot_nsn_vs_z('summary_XKDZ7SI/data/6596_normal/final_state_sim_history.npy', plot_labels=True,  title='$N_{SN} vs z$ [normal SN]', savefig='nsn_z_normal_labels.png')
    #    plot_nsn_vs_z('summary_XKDZ7SI/data/6595_faint/final_state_sim_history.npy',  plot_labels=False, title='$N_{SN} vs z$ [faint SN]', savefig='nsn_z_faint_nolabels.png')
    #    plot_nsn_vs_z('summary_XKDZ7SI/data/6595_faint/final_state_sim_history.npy',  plot_labels=True,  title='$N_{SN} vs z$ [faint SN]', savefig='nsn_z_faint_labels.png')
    
