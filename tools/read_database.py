#!/usr/bin/env python

import os
import os.path as op

import logging
logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',
                    level=logging.INFO)
import argparse
import sqlite3
import numpy as np
import croaks 


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


class WPC_FORMAT:
    table_name = 'SummaryAllProps'
    keymap = {'observationStartMJD': 'mjd',
              'filter': 'band',
              'visitExposureTime': 'exptime',
              'skyBrightness': 'sky',
              'fieldRA': 'Ra',
              'fieldDec': 'Dec',}

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
        
    def get_data(self, cols=None, sql=None, to_degrees=False, new_col_names=None):
        """
        Get the contents of the SQL database dumped into a numpy rec array
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

        logging.info('converting dump into numpy array')
        colnames = [str(c) for c in cols]
        d = np.rec.fromrecords(rows, names=colnames)
        logging.info('update kAtm values')
        katm = np.zeros(len(d), dtype='float64')        
        logging.info('extending output array')
        d = croaks.rec_append_fields(d, data=[katm], names=['kAtm'])
        for b in kAtm.keys():
            d['kAtm'][d['filter'] == b] = kAtm[b]
        logging.info('done.')

        if to_degrees:
            d['fieldRA'] *= (180. / np.pi)
            d['fieldDec'] *= (180. / np.pi)

        if new_col_names is not None:
            self.update_col_names(d, new_col_names)
            
        return d

    def update_col_names(self, d, new_col_names):
        names = list(d.dtype.names)
        d.dtype.names = [new_col_names[n] if n in new_col_names else n for n in d.dtype.names]
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
    parser.add_argument('sql_database',
                        help='cadence database')
    parser.add_argument('--to_degrees', default=False,
                        help='convert Ra and Dec to degrees')
    parser.add_argument('--wpc', default=False,
                        help='.db file uses the format of the white paper call')
    parser.add_argument('--slair', default=False,
                        help='.db file uses the format of the SLAIR (feature based) cadences')
    
    args = parser.parse_args()
    
    print args
    if args.wpc:
        format = WPC_FORMAT
    elif args.slair:
        format = SLAIR_FORMAT
    else:
        format = DEFAULT_FORMAT
        
    reader = Read_Sqlite(args.sql_database, observation_table_name=format.table_name)
    sql = reader.sql_selection(**args.__dict__)
    data = reader.get_data(cols=None, sql=sql,
                           to_degrees=args.to_degrees,
                           new_col_names=format.keymap)
    np.save(args.output, data)


    
