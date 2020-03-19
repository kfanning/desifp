# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 15:14:18 2020

@author: Duan Yutong (dyt@physics.bu.edu)
"""

from bokeh.layouts import row
from bokeh.plotting import figure
from bokeh.models import (
    LinearColorMapper, ColorBar, AdaptiveTicker, LabelSet)
from bokeh.models import ColumnDataSource
from bokeh.models.widgets.tables import (
    DataTable, TableColumn, IntEditor
)
from bokeh.io import curdoc


import os
import pickle
from glob import glob
import posconstants as pc

class PosCalManager:

    table_filename = 'poscal_index.csv'

    def __init__(self):
        temp = os.path.join(pc.dirs['kpno'], '*/*/calibrationdata.pkl')
        self.data_paths = [p for p in glob(temp) if
                           'arc_calibration' in p or 'grid_calibration' in p]

    def init_table(self):
        # if table already exists, load table
        if os.path.isfile(self.table_filename):
            self.df = pd.read_csv(self.table_filename)
            self.table = df.to_dict()  # dict with column names as keys
            self.data_paths = [path for path in self.data_paths
                               if path not in self.table['path']]
        else:  # construct new table
            self.table = d = {'UTC': [], 'expid': [],
                              'test name': [], 'dome': [], 'path': []}
        for path in self.data_paths:  # fill in records for all paths
            with open(os.path.join(path), 'rb') as h:
                data = pickle.load(h)
            d['UTC'].append(data.t_i.isoformat())
            d['expid'].append(data.expid)
            d['test name'].append(data.test_name)
            d['path'].append(path)
        self.df = pd.DataFrame(self.table)
        self.df.to_csv(self.table_filename)


if __name__ == '__main__':
    pcm = PosCalManager()
    print(pcm.df)

