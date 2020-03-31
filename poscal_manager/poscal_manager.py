# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 16:59:03 2020

@author: Duan Yutong (dyt@physics.bu.edu)


Description of environmental conditions that are recorded in table

* Primary Mirror Thermal Control (PMTC):
    tl;dr:
        * only during the day with dome closed when telescope is not in use
        * mirror covers should be closed but could be open
    This draws dome air in through a filter (you can see the large hose on
    the SE side just under the railing on the C floor) into a heat
    exchanger located behind the south journal which valves building
    glycol to cool the air depending on the primary temperature setting
    and dew point.  The air is then blown up through the East elliptical
    tube into the space between the primary mirror and mirror covers.
    This is used only during the day when the telescope is not being used.
    One can hear the vibration from the large blower when it is on.

* Primary Cell Vent Fans (PCVF):
    tl;dr:
        * night time dome closed only
        * mirror covers should be open
    These are the four squirrel cage fans on the center box section of the
    telescope at the NE, SE, SW, and NW locations.  These are controlled
    by the OA and are used to pull air from the space in front of the
    primary around the periphery of the primary mirror and out into the
    dome. The purpose of these is to avoid stagnant air cells which can
    form in the closed space within the center section in front of the
    primary and improve the image quality. As noted below, imaging
    experiments have not suggested any degradation of the image quality
    (and things should actually be better), but as far as I know, no high
    frequency measurements have ever been carried out. Whether these can
    induce vibration of the FVC is an open question. As noted, the OA can
    turn these on or off as desired.  They are used at night with the
    mirror covers open; the PMTC is always turned off when observing.

* Dome Louvers:
    While not properly fans, these are opened to flush air
    through the dome to minimize temperature differences between the telescope
    structure and the ambient air.  Depending on the wind speed and direction,
    it may be possible to set up resonances which could vibrate the telescope
    structure; this could be tested by opening/closing the louvers.

* Dome Fans:
    There are two large fans on the SSW side of the dome on the
    M floor level.  Prior to the installation of the louvers, this provided
    the only method of flushing air through the dome, and they were effectively
    useless in that respect (we carried out a test with neutral buoyancy
    balloons which demonstrated that turning these noisy fans on and off made
    no difference to the air flow in the dome).  I don't know if they are used
    anymore, but they do generate a bit of noise which could couple into
    telescope vibrations.

* B29 Fan:
    This may or may not be re-installed. This is a huge fan which
    sits on the M floor to the East of the hatch and, as its name suggests, is
    extremely noisy. It generates a vertical flow within the closed dome
    during the day to avoid thermal stratification of the air and is supposed
    to facilitate flushing of the dome air once the dome and louvers are opened
    for nighttime observing. The noise makes any work within the dome during
    the day impossible, and this has typically been turned off when any daytime
    work is going on or the mirror covers are open (because of dust that could
    be stirred up).

"""

import sys
import os
import pickle
from glob import glob
from datetime import timezone
import numpy as np
import pandas as pd
from DOSlib.positioner_index import PositionerIndex
sys.path.append('/data/focalplane/pecs/')
if os.environ.get('POSITIONER_LOGS_PATH') is None:
    os.environ['POSITIONER_LOGS_PATH'] = '/data/focalplane/logs'
if os.environ.get('FP_SETTINGS_PATH') is None:
    os.environ['FP_SETTINGS_PATH'] = '/data/focalplane/fp_settings'
import posconstants as pc
from petaltransforms import PetalTransforms


def get_positioner_index():  # get obsXY data for each device
    pi = PositionerIndex()
    pi_df = pd.DataFrame(pi.data).set_index('DEVICE_ID')
    pi_df.columns = pi_df.columns.str.lower()
    pi_df.insert(0, 'device_id', pi_df.index)
    pi_df['obs_x'], pi_df['obs_y'] = np.nan, np.nan
    ptlXYZ_df = pd.read_csv(os.path.join(
            os.getenv('PLATE_CONTROL_DIR',
                      '/software/products/plate_control-trunk'),
            'petal', 'positioner_locations_0530v14.csv'),
        usecols=['device_location_id', 'X', 'Y', 'Z'],
        index_col='device_location_id')
    ptlXYZ_df.index.rename('device_loc', inplace=True)
    for pcid in range(10):
        trans = PetalTransforms(gamma=np.pi/5*(pcid-3))
        obsXY = trans.ptlXYZ_to_obsXYZ(ptlXYZ_df.T.values)[:2, :]
        xy_df = pd.DataFrame(data=obsXY.T, index=ptlXYZ_df.index,
                             columns=['obs_x', 'obs_y'])
        xy_df['device_id'] = (pi_df[pi_df['petal_loc'] == pcid]
                              .set_index('device_loc')['device_id'])
        xy_df.set_index('device_id', inplace=True)
        xy_df = xy_df[xy_df.index.notnull()]
        pi_df.loc[xy_df.index, ['obs_x', 'obs_y']] = xy_df
    return pi_df


class PosCalManager:

    table_path = os.path.join(os.path.abspath('.'), 'poscal_manager',
                              'poscal_index.csv')

    def __init__(self):
        self.pcids = list(range(10))  # initially checked pcids
        temp = os.path.join(pc.dirs['kpno'], '*/*/calibrationdata.pkl')
        self.data_paths = [p for p in glob(temp) if
                           'arc_calibration' in p or 'grid_calibration' in p]
        self._data = {}
        self.init_table()

    def init_table(self):
        cols = ['dome', 'zenith angle', 'tracking', 'PMTC', 'PCVF',
                'dome louvers', 'dome fans', 'B29 fan',
                'cage baffle', 'barrel insulation', 'operator']
        if os.path.isfile(self.table_path):  # csv already exists, load df
            self.df = pd.read_csv(self.table_path, index_col=0)
            self.df_to_table_dict()
            # filter out expids that have are already in the table
            self.data_paths = [path for path in self.data_paths
                               if path not in self.table['data path']]
        else:  # no csv file exists, construct new table dict
            self.table = d = {'UTC': [], 'expid': [], 'test name': [],
                              'data path': []}
            d.update({col: [] for col in cols})
        if self.data_paths:  # fill table dict
            for path in self.data_paths:  # fill in records for all paths
                self.read_pickle(path)
                data = self._data[path]
                if data.posids[0][0] == 'D':
                    continue  # filter out sim runs with Dxxxxx positioners
                d['UTC'].append(data.t_i.astimezone(timezone.utc).strftime(
                    pc.timestamp_format))
                d['expid'].append(data.expid)
                d['test name'].append(data.test_name)
                d['data path'].append(path)
                for col in cols:
                    d[col].append('?')
        self.table_dict_to_df()
        self.len = len(self.df)
        # select the latest arc calibration so that all plots are populated
        self.i_selected = (self.df[self.df['test name'].astype(str)
                           .str.contains('arc')].index[-1])

    def df_to_table_dict(self):
        self.table = self.df.to_dict(orient='list')

    def table_dict_to_df(self):
        self.df = pd.DataFrame(self.table)
        self.df.to_csv(self.table_path)

    def read_pickle(self, path):
        if path not in self._data:
            with open(os.path.join(path), 'rb') as h:
                self._data[path] = pickle.load(h)
        return self._data[path]


if __name__ == '__main__':
    pcm = PosCalManager()
    print(pcm.df)
