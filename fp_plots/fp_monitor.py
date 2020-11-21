# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 12:14:24 2020

@author: Duan Yutong (dyt@physics.bu.edu)

the Singleton is impelemted for DB connection only and FPMonitor re-uses
existing connection. The singleton cannot be applied to FPMonitor class
as Bokeh requires document lock for each document/socket session, requiring
an independent FPMonitor class instance
"""

from collections import OrderedDict
import os
import logging
from datetime import datetime, timedelta
import numpy as np
import pandas as pd
from DOSlib.positioner_index import PositionerIndex
import posconstants as pc
from petaltransforms import PetalTransforms
from db_connect import DBConn
from bokeh.models import ColumnDataSource


class FPMonitor:

    ptlids = {0: '04', 1: '05', 2: '06', 3: '03', 4: '08',
              5: '10', 6: '11', 7: '02', 8: '07', 9: '09'}
    spectro_sns = {0: 'SM04', 1: 'SM10', 2: 'SM05', 3: 'SM06', 4: 'SM01',
                   5: 'SM09', 6: 'SM07', 7: 'SM08', 8: 'SM02', 9: 'SM03'}
    spectro_lns = {0: 'SP0', 1: 'SP1', 2: 'SP2', 3: 'SP3', 4: 'SP4',
                   5: 'SP5', 6: 'SP6', 7: 'SP7', 8: 'SP8', 9: 'SP9'}
    petal_locs = range(10)
    status_colors = OrderedDict([
        ('pospwr_ps1_en', {0: 'crimson', 1: 'limegreen'}),
        ('pospwr_ps2_en', {0: 'crimson', 1: 'limegreen'}),
        ('buff_en1', {0: 'crimson', 1: 'limegreen'}),
        ('buff_en2', {0: 'crimson', 1: 'limegreen'}),
        ('gfapwr_en', {0: 'crimson', 1: 'limegreen'}),
        ('ccdbiasenabled', {0: 'crimson', 1: 'limegreen'}),
        ('gfa_fan_in_en', {0: 'crimson', 1: 'limegreen'}),
        ('gfa_fan_out_en', {0: 'crimson', 1: 'limegreen'}),
        ('fxc_okay', {0: 'red', 1: 'limegreen'}),
        ('pospwr_ps1_fbk', {0: 'crimson', 1: 'limegreen'}),
        ('pospwr_ps2_fbk', {0: 'crimson', 1: 'limegreen'}),
        ('telemetryfault', {0: 'limegreen', 1: 'yellow'}),
        ('canbusfault', {0: 'limegreen', 1: 'yellow'})])

    def __init__(self, refresh_interval=10):
        self.logger = logging.getLogger('FP Telemetry Monitor')
        self.logger.setLevel(logging.WARN)  # log everything, DEBUG level up
        logpath = os.environ.get('DOS_FPPLOTS_LOGS', os.path.abspath('.'))
        log_path = os.path.join(
            logpath, 'fp_plots.log')
        open(log_path, 'a').close()
        fh = logging.FileHandler(log_path, mode='a', encoding='utf-8')
        fh.setFormatter(logging.Formatter(  # log format for each line
            fmt='%(asctime)s %(name)s [%(levelname)-8s]: %(message)s',
            datefmt=pc.timestamp_format))
        self.logger.addHandler(fh)
        self.logger.info('New FPMonitor intance...')
        self.pi = PositionerIndex()
        pi_df = pd.DataFrame(self.pi.data).set_index('DEVICE_ID')
        pi_df.columns = pi_df.columns.str.lower()
        pi_df.insert(0, 'device_id', pi_df.index)
        cols = ['spectro_sn', 'spectro_ln',
                'time_recorded', 'temp', 'posfid_state', 'obs_x', 'obs_y',
                'temp_color', 'line_color']
        dtypes = [str, str, 'datetime64[ns]', np.float32, str,
                  np.float32, np.float32, np.float32, np.float32, str]
        data = {col: pd.Series(dtype=dt) for col, dt in zip(cols, dtypes)}
        # import pdb; pdb.set_trace()
        self.data = pi_df.join(pd.DataFrame(data=data))
        # add dummy centre cap ring so data won't be all NaN or bokeh crashes
        # self.data.loc['centre cap ring', ['device_loc', 'device_type', 'temp_color', 'line_color']] = -1, 'centre cap ring', 0, 'white'
        # self.data.loc['centre cap ring', ''] = 0
        # self.data.loc['centre cap ring', 'line_color'] = 'white'
        # create masks for pos and fid
        self.posmask = ((self.data['device_type'] == 'POS') |
                        (self.data['device_type'] == 'ETC'))
        self.fidmask = ((self.data['device_type'] == 'GIF')
                        | (self.data['device_type'] == 'FIF'))
        self.data.loc[self.posmask, 'line_color'] = 'white'
        self.data.loc[self.fidmask, 'line_color'] = 'green'
        # add obsXYZ positions to data for plotting
        path = os.path.join(
            os.getenv('PLATE_CONTROL_DIR',
                      '/software/products/plate_control-trunk'),
            'petal', 'positioner_locations_0530v18.csv')
        ptlXYZ_df = pd.read_csv(path,
                                usecols=['device_location_id', 'X', 'Y', 'Z'],
                                index_col='device_location_id')
        ptlXYZ_df.index.rename('device_loc', inplace=True)
        self.data_status = pd.DataFrame(
            data={'x': 21.5 + 9.5 * (
                        np.arange(len(self.status_colors.keys()))//9),
                  'y': np.power(10, 2-0.17*(
                        np.arange(len(self.status_colors.keys())) % 9))},
            index=FPMonitor.status_colors.keys())
        self.data_status.index.rename('status_field', inplace=True)
        for petal_loc in self.petal_locs:
            trans = PetalTransforms(gamma=np.pi/5*(petal_loc-3))
            obsXY = trans.ptlXYZ_to_obsXYZ(ptlXYZ_df.T.values)[:2, :]
            xy_df = pd.DataFrame(data=obsXY.T, index=ptlXYZ_df.index,
                                 columns=['obs_x', 'obs_y'])
            xy_df['device_id'] = (
                self.data[self.data['petal_loc'] == petal_loc]
                .set_index('device_loc')['device_id'])
            xy_df.set_index('device_id', inplace=True)
            xy_df = xy_df[xy_df.index.notnull()]
            self.data.loc[xy_df.index, ['obs_x', 'obs_y']] = xy_df
            # add spectro serial and logical numbers for each petal location
            self.data.loc[self.data['petal_loc'] == petal_loc,
                          'spectro_sn'] = self.spectro_sns[petal_loc]
            self.data.loc[self.data['petal_loc'] == petal_loc,
                          'spectro_ln'] = self.spectro_lns[petal_loc]
            # pc telemetry status data
            self.data_status[f'val_{petal_loc}'] = np.nan
            self.data_status[f'color_{petal_loc}'] = 'grey'
        self.source = ColumnDataSource(self.data)
        self.source_status = ColumnDataSource(self.data_status)
        # self.delta_t = timedelta(seconds=refresh_interval)
        # self.last_updated = pc.now() - self.delta_t
        # load past 5 min and use the latest readout
        self.update_data_and_sources(update_hist=False)

    def query_db(self, minutes=5, table='pc_telemetry_can_all'):
        time_range = timedelta(minutes=minutes)
        time_cut = (datetime.utcnow() - time_range).strftime(
            '%Y-%m-%d %H:%M:%S')
        query = pd.read_sql_query(
            f"""SELECT * FROM {table}
                WHERE time_recorded >= '{time_cut}'""",
            DBConn.conn).sort_values('time_recorded')  # ascending, old to new
        return query  # return query even if empty and no data in 5 min

    def process_pc_telemetry_can_all(self, query, latest_row_only=True):
        '''adds the result of a sql query to the data source
           11/17/2020 now expects querry of PTL-TEMPS table which contains
                the required data.'''
        for petal_loc in self.petal_locs:
            query_petal = query[query['pcid'] == petal_loc]
            if query_petal.empty:  # no telemetry in the past x min for any ptl
                self.logger.info(
                    f'No temperature telemetry found for PC0{petal_loc}, resetting...')
                mask = self.data['petal_loc'] == petal_loc
                self.data.loc[mask, 'time_recorded'] = pc.timestamp_str()
                self.data.loc[mask, 'temp'] = np.nan
                self.data.loc[mask, 'posfid_state'] = np.nan
            else:
                for i, series in query_petal.iterrows():
                    if latest_row_only and i < len(query_petal)-1:
                        # current row i is not the last/latest row, skip it
                        continue
                    device_ids, temps = [], []
                    # gather all device_id and temp in this row into 2 lists
                    for device_loc, temp in series['posfid_temps'].items():
                        if 'can' in device_loc:
                            self.logger.warning(
                                f"Skipping DB entry in old telemetry format "
                                f"submitted by PC-{series['pcid']} "
                                f"at {series['time']}")
                            break
                        device_id = self.pi.find_by_petal_loc_device_loc(
                            series['pcid'], device_loc, key='DEVICE_ID')
                        device_ids.append(device_id)
                        temps.append(temp)
                    self.logger.debug(f'Processed {len(device_ids)} devices '
                                      f'for PC0{series.pcid}')
                    #time = series['time_recorded'].strftime(
                    #    '%Y-%m-%dT%H:%M:%S%z')
                    dtime = series['time']
                    if isinstance(dtime, pd.Timestamp):
                        dtime = dtime.isoformat()
                    # New format returns isoformat string. Need to cut around second decimals
                    time = dtime[:19] + dtime[-6:]
                    self.data.loc[device_ids,
                                  'time_recorded'] = time #time[:-2]+':'+time[-2:]
                    self.data.loc[device_ids, 'temp'] = temps
                    self.data.loc[device_ids,
                                  'posfid_state'] = series['posfid_state']
        # update temperature for each device_type: color = temp - temp_mean
        # for mask in [self.posmask, self.fidmask]:
        #     self.data.loc[mask, 'temp_color'] = (
        #         self.data[mask]['temp'] - self.data[mask]['temp'].mean())
        self.data['temp_color'] = self.data['temp']

    def process_pc_telemetry_status(self, query):
        '''
        Reads telemetry_status table from PC telemetry
        11/17/2020 changed to PTL-STATUS table which is more or less equivalent
        '''
        for petal_loc in self.petal_locs:
            query_petal = query[query['pcid'] == petal_loc]
            if query_petal.empty:
                self.logger.info(
                    f'No PTL-STATUS for PC0{petal_loc}, resetting...')
                for field in self.status_colors.keys():
                    self.data_status.loc[field, f'val_{petal_loc}'] = np.nan
                    self.data_status.loc[field, f'color_{petal_loc}'] = 'grey'
            else:
                series = query_petal.iloc[-1]  # the last row, latest entry
                for field in self.status_colors.keys():
                    value = series[field]
                    # cast values to binary int
                    if value == 0 or value == 1:
                        pass
                    elif value is None or False:
                        value = 0
                    elif value is True:
                        value = 1
                    elif np.isnan(value):
                        value = 0
                    else:
                        self.logger.error(
                            f'Bad status value: {field} = {value}, '
                            f'type {type(value)}')
                        continue  # skip bad vaue for this field
                    self.data_status.loc[field,
                                         f'val_{petal_loc}'] = series[field]
                    self.data_status.loc[field, f'color_{petal_loc}'] = (
                            self.status_colors[field][value])

    def update_data_and_sources(self, minutes=5, update_hist=True):
        self.logger.info('Data source update requested')
        # if pc.now() - self.last_updated < self.delta_t:
        #     self.logger.info('update too frequent, skipping...')
        #     return  # too frequent, wait until delta_t has elapsed
        self.logger.info(f'Updating plotting sources with new telemetry data '
                         f'in the past {minutes} minutes...')
        # do one SQL query for all telemetry data in the past x min
        query = self.query_db(minutes=minutes, table='pc_ptl_temps')
        self.process_pc_telemetry_can_all(query)
        query = self.query_db(minutes=minutes, table='pc_ptl_status')
        self.process_pc_telemetry_status(query)
        self.source.data = self.data
        self.source_status.data = self.data_status
        self.generate_hist_data()
        if hasattr(self, 'source_hist'):
            self.source_hist.data = self.data_hist
        else:
            self.source_hist = ColumnDataSource(self.data_hist)
        # self.last_updated = pc.now()
        return True

    def generate_hist_data(self, n_bins=10):
        '''requires updated data; creates source for histogram, fields are
        left_0_pos, left_0_fid, ..., right_0_pos, right_0_fid, ...,
        top_0_pos,  bottom_0_pos,...'''
        self.data_hist = {}
        for petal_loc in self.petal_locs:
            petal_data = self.data[self.data['petal_loc'] == petal_loc]
            posmask = petal_data['device_type'] == 'POS'
            for device_type, mask in zip(['pos', 'fid'], [posmask, ~posmask]):
                series = petal_data[mask]['temp']
                hist, edges = np.histogram(series[series.notnull()],
                                           density=False, bins=n_bins)
                if np.all(hist == 0):  # empty telemetry, hist is all zeros
                    # use alternative top so bokeh histogram doesn't crash
                    hist = hist.astype(np.float64) + 0.1
                self.data_hist[f'top_{petal_loc}_{device_type}'] = hist
                self.data_hist[f'left_{petal_loc}_{device_type}'] = edges[:-1]
                self.data_hist[f'right_{petal_loc}_{device_type}'] = edges[1:]
