#!/usr/bin/env python

import os
import os.path as op
import sys
import numpy as np
import pylab as pl
import re
import matplotlib.cm as cm
import glob
import json


def nsn_and_median_zmax(dirname):
    zmax = np.load(dirname + os.sep + 'zmax_inst_history.npy')
    nsn = np.load(dirname + os.sep + 'nsn_tot_history.npy')
    return np.median(zmax['val']), nsn['val'].max()


pattern = re.compile('(\d+)_(.+)_(\d{5})_0_(\d{5})_0_(\d+)_(faint|normal)')

family_color = {
    '2018-01': cm.jet(10),
    'WPC': cm.jet(55),
    'SLAIR': cm.jet(210),
    'altsched': cm.jet(240),
    }


class PlotData():
    
    def __init__(self, filename):
        """
        collect the result directories, and the plot options.
        """
        self.root_dir = op.dirname(filename)
        with open(filename) as f:
            d = json.load(f)
            self.options = d['cadences']
            self.default = d['default']
            
        rr = []
        for dn in glob.glob(self.root_dir + os.sep + '*faint') + \
                  glob.glob(self.root_dir + os.sep + '*normal'):
            r = pattern.search(dn)
            if r is None:
                logging.warning('unable to parse: %s' % dn)
                continue
            _, cadence, _, _, nside, sntype = r.groups()
            opts = self.options.get(cadence, self.default)
            rr.append((cadence, nside, sntype, dn, opts['draw'], opts['family'], opts['label'], opts['legend'], opts['rolling'], opts['marker'], opts['annotate'], opts['name']))
            
        self.db = np.rec.fromrecords(rr, names=['cadence', 'nside', 'sntype', 'dirname', 'draw', 'family', 'label', 'legend', 'rolling', 'marker', 'annotate', 'name'])
        self.collect()

    def collect(self):
        """
        """
        nsn_zmax = [nsn_and_median_zmax(d.dirname) for d in self.db]
        self.nsn_zmax = np.rec.fromrecords(nsn_zmax, names=['zmax', 'nsn'])
        
    def plot(self, sntype='normal', xlim=(0.25, 0.55), ylim=(0., 400000.), xlabel='$z_{med}$', ylabel='$N_{z<z_{med}}$'):
        """
        """
        fig = pl.figure(figsize=(16,9))
        idx = self.db['sntype'] == sntype
        data = self.nsn_zmax[idx]
        db = self.db[idx]
        pl.plot(data['zmax'], data['nsn'],
                markeredgecolor='w', marker='.', markerfacecolor='w',
                ls='')
        
        dx, dy = 0.005, -3000
        for d,o in zip(data,db):
            if not o['draw']:
                continue
            label = o['family'] if o['label'] else None
            if o['rolling']:
                pl.plot(d['zmax'], d['nsn'],
                        color=family_color[o['family']], markerfacecolor='w', marker=o['marker'], markersize=11, # marker=o['marker'], 
                        ls='', label=label)            
            else:
                pl.plot(d['zmax'], d['nsn'],
                        color=family_color[o['family']], marker=o['marker'], markersize=11, # marker=o['marker']
                        ls='', label=label)
            if o['annotate']:
                pl.text(d['zmax']+dx, d['nsn']+dy, o['name'], color=family_color[o['family']], fontsize=11)

        ax = pl.gca()
        ax.tick_params(axis='both', labelsize=14)
        pl.xlabel(xlabel, fontsize=18)
        pl.ylabel(ylabel, fontsize=18)
        pl.xlim(xlim)
        pl.ylim(ylim)
        pl.legend(loc='upper left', ncol=2, shadow=True)
        pl.grid(ls=':')
        pl.subplots_adjust(left=0.15, right=0.95, top=0.95)

        
        

#if __name__ == "__main__":
if False:
    plot_opts = PlotOptions('./plot_options.json')
    dx, dy = 0.005, -3000
    
    data = []
    for nm in sys.argv[1:]:
        r = pattern.search(nm)
        _, cadence, mjd_min, mjd_max, _, which_sn = r.groups()
        cad, col, marker, plot_cadence, plot_label, annotate, rolling, family = plot_opts(cadence)
        try:
            data.append(nsn_and_median_zmax(nm) + (cadence, cad, col, marker, family, plot_cadence, plot_label, annotate, rolling))            
        except:
            print ' (*) unable to process: ', nm
            print sys.exc_info()
            continue
    data = np.rec.fromrecords(data, names=['zmax', 'nsn', 'cadence', 'cad', 'color', 'marker', 'family', 'plot_cadence', 'plot_label', 'annotate', 'rolling'])
    print data 

    fig = pl.figure(figsize=(16,9))
    pl.plot(data['zmax'], data['nsn'],
                markeredgecolor='w', marker=',', markerfacecolor='w',
                 ls='')
    for d in data:
        if not d['plot_cadence']:
            continue
        label = d['family'] if d['plot_label'] else None
        if d['rolling']:
            pl.plot(d['zmax'], d['nsn'],
                    color=family_color[d['family']], markerfacecolor='w', marker=d['marker'], markersize=11,
                    ls='', label=label)            
        else:
            pl.plot(d['zmax'], d['nsn'],
                    color=family_color[d['family']], marker=d['marker'], markersize=11,
                    ls='', label=label)
        if d['annotate']:
            pl.text(d['zmax']+dx, d['nsn']+dy, d['cad'], color=family_color[d['family']], fontsize=14)
        
    ax = pl.gca()
    ax.tick_params(axis='both', labelsize=14)
    pl.xlabel('$z_{max}$', fontsize=18)
    pl.ylabel('$Number\ of\ well sampled\ SNe$', fontsize=18)
    pl.xlim((0.25, 0.55))
    pl.ylim((0., 400000.))
    pl.legend(loc='upper left', ncol=2, shadow=True)
    pl.grid(ls=':')
    pl.subplots_adjust(left=0.15, right=0.95, top=0.95)
    
    fig.savefig('cadence_summary_plot.png', bbox_inches='tight')










# class PlotOptions:

#     default_plot_opts = {
#         'feature_rolling_half_mask_10yrs': ('Feature Rolling 1/2', 'y', 'o',      1,  0,  1, True, '2018-01'),
#         'minion_1016': ('Minion 1016',                             'y', 'o',      1,  1,  1, False, '2018-01'),
#         'feature_rolling_twoThird_10yrs': ('Feature Rolling 2/3',  'y', 'o',      1,  0,  1, True,  '2018-01'),
#         'feature_baseline_10yrs': ('Feature Baseline',             'y', 'o',      1,  0,  1, False,  '2018-01'),
        
#         'alt_sched': ('AltSched',                                  'b',    's',   0,  1,  1, False, 'altsched'),
#         'alt_sched_rolling': ('AltSched Rolling',                  'b',    's',   0,  0,  1, True, 'altsched'),
#         'altsched_good_weather': ('AltSched',                      'b',    's',   1,  1,  1, False, 'altsched'),
#         'altsched_rolling_good_weather': ('AltSched Rolling',      'b',    's',   1,  0,  1, True, 'altsched'),
#         'altsched_18__90_30': ('AltSched wide',                    'b',    's',   1,  0,  1, False, 'altsched'),
#         'altsched_18__90_40': ('AltSched wide',                    'b',    's',   1,  0,  0, False, 'altsched'),
        
#         'mothra_2045':   ('Mothra 2045',    'gray',  'v',   1,   0,  1, True,  'WPC'), # rolling
#         'kraken_2035':   ('Kraken 2035',    'gray',  'v',   1,   1,  0, False, 'WPC'),
#         'colossus_2667': ('Colossus 2667',  'gray',  'v',   1,   0,  1, False, 'WPC'),
#         'colossus_2664': ('Colossus_ 2664', 'gray',  'v',   1,   0,  0, False, 'WPC'),
#         'colossus_2665': ('Colossus 2665',  'gray',  'v',   1,   0,  0, False, 'WPC'),
#         'pontus_2002':   ('Pontus 2002',    'gray',  'v',   1,   0,  1, False, 'WPC'),
#         'kraken_2026':   ('Kraken 2026',    'gray',  'v',   1,   0,  0, False, 'WPC'),
#         'baseline2018a': ('Baseline 2018a', 'gray',  'v',   1,   0,  0, False, 'WPC'),
#         'pontus_2489':   ('Pontus 2489',    'gray',  'v',   1,   0,  1, False, 'WPC'),
#         'pontus_2502':   ('Pontus 2502',    'gray',  'v',   1,   0,  1, True,  'WPC'), # rolling
#         'kraken_2036':   ('Kraken 2036',    'gray',  'v',   1,   0,  1, True,  'WPC'), # rolling

#         'kraken_2042':   ('Kraken 2042',    'red',   'o',   1,   0,  1, False,  'WPC'), # 
#         'kraken_2044':   ('Kraken 2044',    'red',   'o',   1,   0,  1, False,  'WPC'), # 
#         'mothra_2049':   ('Mothra 2049',    'red',   'o',   1,   0,  1, True,  'WPC'), # rolling
#         'nexus_2097':    ('Nexus 2097',     'red',   'o',   1,   0,  1, True,  'WPC'), # rolling

#         'blobs_mix_zmask10yrs':  ('Slair mix zmask',  'g', '>',         0, 0, 1, False, 'SLAIR'),
#         'blobs_same_10yrs':      ('Slair same',       'g', '>',         0, 0, 1, False, 'SLAIR'),
#         'blobs_same_zmask10yrs': ('Slair same zmask',  'g', '>',        0, 0, 1, False, 'SLAIR'),        
#         'rolling_10yrs':         ('Slair rolling',          'g', '>',   1, 0, 1, True,  'SLAIR'), # rolling
#         'rolling_mix_75_10yrs':   ('Slair rolling 75', 'g', '>',        1, 0, 1, True, 'SLAIR'),
#         'cadence_mix_10yrs': ('Slair mix',       'g', '>',              1, 1, 1, False, 'SLAIR'),
#         'rolling_mix_10yrs': ('Slair rolling mix', 'g', '>',            1, 0, 1, True, 'SLAIR'),
#     }

#     def __init__(self):
#         self.plot_opts = self.default_plot_opts

#     def __call__(self, k):
#         return self.plot_opts.get(k, ('?', 'gray', '.', 0, 0, 0, True, '?'))
    
