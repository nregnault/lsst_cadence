#!/usr/bin/env python

import os
import os.path as op
from exceptions import ValueError

import logging
logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',
                    level=logging.INFO)
import argparse
import sqlite3
import numpy as np
import croaks 
from numpy.lib.recfunctions import join_by


kAtm = { 'u': 0.50,
         'g': 0.21,
         'r': 0.13,
         'i': 0.10,
         'z': 0.07,
         'y': 0.18 }


dtype = ['index', 'observationId', 'night', 
         'observationStartTime', 
         'observationStartMJD', 
         'observationStartLST', 
         'filter', 'proposalId', 
         'fieldRA', 'fieldDec', 'angle', 
         'altitude', 'azimuth', 'numExposures', 
         'visitTime', 'visitExposureTime',
         'airmass', 'skyBrightness', 'cloud', 'seeingFwhm500', 
         'seeingFwhmGeom', 'seeingFwhmEff', 'fiveSigmaDepth', 
         'moonRA', 'moonDec', 'moonAlt', 'moonAz', 'moonDistance', 'moonPhase', 
         'sunRA', 'sunDec', 'sunAlt', 'sunAz', 'solarElong', 
         'slewTime', 'slewDistance', 'paraAngle', 'rotTelPos', 'rotSkyPos', ]


# these are the default columns Philippe used to retrieve
# now, by default, we get all the columns
default_cols = ['observationStartMJD', 'fiveSigmaDepth', 'filter', 'observationStartTime',
                'fieldRA', 'fieldDec', 'visitTime', 'visitExposureTime', 'numExposures',
                'seeingFwhm500', 'seeingFwhmGeom', 'seeingFwhmEff', 'moonPhase', 'airmass',
                'skyBrightness']            


def _katm_column(data):
    """
    """
    logging.info('update kAtm values')
    katm = np.zeros(len(data), dtype='float64')        
    for b in kAtm.keys():
        katm[data['filter'] == b] = kAtm[b]
    return katm, 'kAtm'


def _remove_duplicate_entries(data):
    """The cadence database dumps may contain duplicate entries, allocated
    simultaneously to two surveys (typically WFD and DDF). We try to
    remove them here.
    """
    logging.info('removing duplicate entries')
    c = np.bincount(data['observationId'])
    dup_ids = np.where(c>1)
    idx = np.in1d(data['observationId'], dup_ids)
    logging.info('%d duplicates detected' % idx.sum())
    oid = data['observationId'][idx]
    d = np.zeros(len(oid))
    d[:-1] = (oid[1:] - oid[:-1]) == 0
    d *= data['observationId'][idx]
    data['observationId'][idx] = d
    r = data[data['observationId']>0]
    logging.info('size: %d -> %d' % (len(data), len(r)))
    return r
    

def _desc_dither_columns(data, dithers):
    logging.info('adding dithers')
    d = join_by('observationId', data, dithers, jointype='inner', 
                defaults={'descDitheredRA': 0., 'descDitheredDec': 0., 'descDitheredRotTelPos': 0.}, 
                usemask=False)
    #    nm = []
    #    for nn in d.dtype.names:
    #        if nn == 'observationId':
    #            nm.append(nn)
    #        else:
    #            nm.append(nn[:-1])
    #    d.dtype.names = nm
    d['rotTelPos'] = d['descDitheredRotTelPos']
    return d
        
        
class WPC_FORMAT:
    table_name = 'SummaryAllProps'
    keymap = {'observationStartMJD': 'mjd',
              'filter': 'band',
              'visitExposureTime': 'exptime',
              'skyBrightness': 'sky',
              'descDitheredRA': 'Ra',
              'descDitheredDec': 'Dec',
              #              'descDitheredRotTelPos': 'rotTelPos',
    }

class SLAIR_FORMAT:
    table_name = 'observations'
    keymap = {'RA': 'Ra',
              'dec': 'Dec',
              'filter': 'band',
              'FWHMeff': 'seeingFwhmEff',
              'skybrightness': 'sky',
              'fivesigmadepth': 'fiveSigmaDepth',}

class DEFAULT_FORMAT:
    table_name = 'SummaryAllProps'
    keymap = {}


class Read_Sqlite:
    
    def __init__(self, dbfile, **sel):
        """
        """
        self.dbfile = dbfile
        conn = sqlite3.connect(dbfile)
        self.cur = conn.cursor()
        # get the table list
        self.cur.execute("SELECT name FROM sqlite_master WHERE type='table';")
        self.tables = self.cur.fetchall()
        self.sql = self.sql_selection(**sel)
        #        self.data = self.groom_data(self.get_data(sql))
        self.observation_table_name = sel.get('observation_table_name', 'SummaryAllProps')
        self.dithers = sel.get('dithers', None)
        
    def sql_selection(self, **sel):
        sql = ''
        if 'year' in sel and sel['year'] is not None:
            y = sel['year']
            sql += 'night > %i and night < %i' % ((y-1)*365.25, y*365.25)
        if 'mjd_min' in sel and sel['mjd_min'] is not None:
            if len(sql) > 0: sql += ' and '
            sql += 'observationStartMJD > %f' % sel['mjd_min']
        if 'mjd_max' in sel and sel['mjd_max'] is not None:
            if len(sql) > 0: sql += ' and '
            sql += 'observationStartMJD < %f' % sel['mjd_max']            
        if 'proposalId' in sel and sel['proposalId'] is not None:
            if len(sql) > 0: sql += ' and '
            sql += 'proposalId=%d' % sel['proposalId']
        return sql
        
    def _get_rows(self, cols=None, sql=None):
        """
        Get the contents of the SQL database as a list of rows
        """
        sql_request = 'SELECT '
        if cols is None:
            sql_request += ' * '            
            self.cur.execute('PRAGMA TABLE_INFO(%s)' % self.observation_table_name)
            r = self.cur.fetchall()
            cols = [c[1] for c in r]
        else:
            sql_request += ','.join(cols)
            
        sql_request += 'FROM %s' % self.observation_table_name
        if sql is not None and len(sql) > 0:
            sql_request += ' WHERE ' + sql
        sql_request += ';'
        
        logging.info('fetching data from db')
        logging.info('request: %s' % sql_request)
        self.cur.execute(sql_request)
        rows = self.cur.fetchall()
        logging.info('done. %d entries' % len(rows))
        return rows, cols
    
    def _to_recarray(self, rows, cols, new_col_names=None):
        """
        convert list of rows into a recarray 
        """
        logging.info('converting dump into numpy array')
        colnames = [str(c) for c in cols]
        d = np.rec.fromrecords(rows, names=colnames)
        
        # extend the array with additional values
        cols = [_katm_column(d)]
        logging.info('extending output array')
        d = croaks.rec_append_fields(d, data=[c[0] for c in cols], names=[c[1] for c in cols])
        logging.info('done.')
        
        # there may be duplicate entries 
        d = _remove_duplicate_entries(d)

        # add the precomputed dithers
        if self.dithers is not None:
            print len(d), len(self.dithers)
            d = _desc_dither_columns(d, self.dithers)            
            
        # remapping column names to pipeline default format
        if new_col_names is not None:
            self._update_col_names(d, new_col_names)

        return d

    def _update_col_names(self, d, new_col_names):
        names = list(d.dtype.names)
        d.dtype.names = [new_col_names[n] if n in new_col_names else n for n in d.dtype.names]
        return d
        
    def get_data(self, cols=None, sql=None, to_degrees=False, new_col_names=None):
        rows, cols = self._get_rows(cols)
        d = self._to_recarray(rows, cols, new_col_names)
        if to_degrees:
            d['Ra'] *= (180. / np.pi)
            d['Dec'] *= (180. / np.pi)
        return d
            


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Dump a cadence database into a numpy npy file.')
    parser.add_argument('-o', '--output',
                        default='cadence.npy',
                        help='output cadence file')
    parser.add_argument('-b', '--band',
                        default=None, type=str,
                        help='limit the output to the specified band')
    parser.add_argument('--mjd_min',
                        default=None, type=float,
                        help='only pointings taken after mjd')
    parser.add_argument('--mjd_max',
                        default=None, type=float,
                        help='only pointings taken before mjd')
    parser.add_argument('--year',
                        default=None, type=float,
                        help='only pointing taken between Jan 1st and Dec 31st of YEAR')
    parser.add_argument('--csv-dithers',
                        default=None, type=str,
                        help='auxiliary CSV file containing the dithers applied to each exposure')
    parser.add_argument('sql_database',
                        help='cadence database')
    parser.add_argument('--to_degrees', default=False,
                        help='convert Ra and Dec to degrees')
    parser.add_argument('--wpc', default=False, action='store_true',
                        help='.db file uses the format of the white paper call')
    parser.add_argument('--slair', default=False, action='store_true',
                        help='.db file uses the format of the SLAIR (feature based) cadences')
    
    args = parser.parse_args()
    
    print args
    if args.wpc:
        format = WPC_FORMAT
    elif args.slair:
        format = SLAIR_FORMAT
    else:
        format = DEFAULT_FORMAT
        

    if args.csv_dithers is not None:
        logging.info('loading dithers from auxiliary csv file: %s' % args.csv_dithers)
        dithers = croaks.NTuple.fromcsv(args.csv_dithers)
    else:
        dithers = None
        
    reader = Read_Sqlite(args.sql_database, observation_table_name=format.table_name, dithers=dithers)
    sql = reader.sql_selection(**args.__dict__)
    data = reader.get_data(cols=None, sql=sql,
                           to_degrees=args.to_degrees,
                           new_col_names=format.keymap)
    
    np.save(args.output, data)

