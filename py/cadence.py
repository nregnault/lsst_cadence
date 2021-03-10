"""
"""

import os
import os.path as op
import sys



class SummaryData():
    def __init__(self, id, summary_data):
        self.id = id
        self.data = summary_data[id]
        self.summary_data = summary_data

    def nsn(self):
        return self.data['nsn_tot']
    def zmax(self):
        return self.data['zmax']
    def cadence(self):
        return self.data['cadence']
        


class Cadence():
    
    def __init__(self, name, version_tag):
        """
        """
        self.name=name
        self.version_tag = version_tag
        self.id = {}
        self.summary_data = {}
        self.movie_file = {}
        
    @property
    def full_name(self):
        return self.summary_data['normal'].data['cadence_name'].decode()
        
    def nsn(self, sntype):
        return self.summary_data[sntype].nsn()

    def zmax(self, sntype):
        return self.summary_data[sntype].zmax()

    def cadence(self):
        return self.summary_data['normal'].cadence()

