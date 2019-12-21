'''
start server with the following command:
    bokeh serve fp_plots --allow-websocket-origin=desi-2.kpno.noao.edu:5006

view at:
    http://desi-2.kpno.noao.edu:5006/main
'''
from collections import OrderedDict
import os
import logging
from datetime import datetime, timedelta
import numpy as np
import psycopg2
import pandas as pd
import posconstants as pc
from petaltransforms import PetalTransforms
from DOSlib.positioner_index import PositionerIndex
from bokeh.io import curdoc, output_file, save
from bokeh.layouts import row, gridplot
from bokeh.palettes import Magma256
from bokeh.plotting import figure
from bokeh.models import (
    ColumnDataSource, LinearColorMapper, ColorBar, AdaptiveTicker, LabelSet)

ptlids = {0: '04', 1: '05', 2: '06', 3: '03', 4: '08',
          5: '10', 6: '11', 7: '02', 8: '07', 9: '09'}
spectro_sns = {0: 'SM04', 1: 'SM10', 2: 'SM05', 3: 'SM06', 4: 'SM01',
               5: 'SM09', 6: 'SM07', 7: 'SM08', 8: 'SM02', 9: 'SM03'}
spectro_lns = {0: 'SP0', 1: 'SP1', 2: 'SP2', 3: 'SP3', 4: 'SP4',
               5: 'SP5', 6: 'SP6', 7: 'SP7', 8: 'SP8', 9: 'SP9'}
timespec = 'seconds'  # for datetime isoformat
refresh_interval = 10  # seconds
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
logger = logging.getLogger('FP Telemetry Monitor')
logger.setLevel(logging.DEBUG)  # log everything, DEBUG level up
log_path = os.path.join(os.path.abspath('.'), 'fp_plots', 'fp_plots.log')
open(log_path, 'a').close()
fh = logging.FileHandler(log_path, mode='a', encoding='utf-8')
fh.setFormatter(logging.Formatter(  # log format for each line
    fmt='%(asctime)s %(name)s [%(levelname)-8s]: %(message)s',
    datefmt=pc.timestamp_format))
logger.addHandler(fh)


def init_data():
    path = os.getenv('DOS_POSITIONERINDEXTABLE')
    if not os.path.isfile(path):
        raise ValueError(f'$DOS_POSITIONERINDEXTABLE = {path}')
    pi_df = pd.read_csv(path)
    pi_df.columns = pi_df.columns.str.lower()
    pi_df.set_index('device_id', inplace=True)
    pi_df.insert(0, 'device_id', pi_df.index)
    cols = ['spectro_sn', 'spectro_ln',
            'time_recorded', 'temp', 'posfid_state', 'obs_x', 'obs_y',
            'temp_color', 'line_color']
    dtypes = [str, str, 'datetime64[ns]', np.float32, str,
              np.float32, np.float32, np.float32, np.float32, str]
    data = {col: pd.Series(dtype=dt) for col, dt in zip(cols, dtypes)}
    data = pi_df.join(pd.DataFrame(data=data, index=pi_df.index))
    # add dummy centre cap ring so data won't be all NaN or bokeh crashes
    data.loc['centre cap ring', 'device_type'] = 'centre cap ring'
    data.loc['centre cap ring', 'temp_color'] = 0
    data.loc['centre cap ring', 'line_color'] = 'white'
    # create masks for pos and fid
    posmask = data['device_type'] == 'POS'
    fidmask = (data['device_type'] == 'GIF') | (data['device_type'] == 'FIF')
    data.loc[posmask, 'line_color'] = 'white'
    data.loc[fidmask, 'line_color'] = 'green'
    # add obsXYZ positions to data for plotting
    path = os.path.join(os.getenv('PLATE_CONTROL_DIR',
                                  '/software/products/plate_control-trunk'),
                        'petal', 'positioner_locations_0530v14.csv')
    ptlXYZ_df = pd.read_csv(path,
                            usecols=['device_location_id', 'X', 'Y', 'Z'],
                            index_col='device_location_id')
    ptlXYZ_df.index.rename('device_loc', inplace=True)
    data_status = pd.DataFrame(
        data={'x': 22 + 9 * (np.arange(len(status_colors.keys()))//9),
              'y': np.power(
                      10, 2-0.17*(np.arange(len(status_colors.keys())) % 9))},
        index=status_colors.keys())
    data_status.index.rename('status_field', inplace=True)
    for petal_loc in petal_locs:
        trans = PetalTransforms(gamma=np.pi/5*(petal_loc-3))
        obsXY = trans.ptlXYZ_to_obsXYZ(ptlXYZ_df.T.values)[:2, :]
        xy_df = pd.DataFrame(data=obsXY.T, index=ptlXYZ_df.index,
                             columns=['obs_x', 'obs_y'])
        # use device_loc for update as some are defined but not installed
        petal_data = (data[data['petal_loc'] == petal_loc]
                      .reset_index().set_index('device_loc'))
        petal_data.update(xy_df)
        data.update(petal_data.reset_index().set_index('device_id'))
        # add spectro serial and logical numbers for each petal location
        data.loc[data['petal_loc'] == petal_loc,
                 'spectro_sn'] = spectro_sns[petal_loc]
        data.loc[data['petal_loc'] == petal_loc,
                 'spectro_ln'] = spectro_lns[petal_loc]
        # pc telemetry status data
        data_status[f'val_{petal_loc}'] = np.nan
        data_status[f'color_{petal_loc}'] = 'grey'
    return (data, ColumnDataSource(data),
            data_status, ColumnDataSource(data_status),
            posmask, fidmask)


def query_db(minutes=5, table='pc_telemetry_can_all'):
    time_range = timedelta(minutes=minutes)
    time_cut = (datetime.utcnow() - time_range).strftime('%Y-%m-%d %H:%M:%S')
    query = pd.read_sql_query(
        f"""SELECT * FROM {table}
            WHERE time_recorded >= '{time_cut}'""",
        conn).sort_values('time_recorded')  # ascending order, old to n
    # return query_db(minutes=minutes+1, table=table) if query.empty else query
    return query  # just return query even if it is empty and no data in 5 min


def process_pc_telemetry_can_all(query, latest_row_only=True):
    '''adds the result of a sql query to the data source'''
    for petal_loc in petal_locs:
        query_petal = query[query['pcid'] == petal_loc]
        if query_petal.empty:  # no telemetry in the past x min for any petal
            logger.info(
                f'No CAN telemetry found for PC0{petal_loc}, resetting...')
            mask = data['petal_loc'] == petal_loc
            data.loc[mask, 'time_recorded'] = pc.timestamp_str()
            data.loc[mask, 'temp'] = np.nan
            data.loc[mask, 'posfid_state'] = np.nan
        else:
            for i, series in query_petal.iterrows():
                if latest_row_only and i < len(query_petal)-1:
                    # current row i is not the last/latest row, skip it
                    continue
                device_ids, temps = [], []
                # gather all device_id and temp info in this row into 2 lists
                for device_loc, temp in series['posfid_temps'].items():
                    if 'can' in device_loc:
                        logger.warning(
                            f"Skipping DB entry in old telemetry format "
                            f"submitted by PC-{series['pcid']} "
                            f"at {series['time_recorded']}")
                        break
                    device_id = pi.find_by_petal_loc_device_loc(
                        series['pcid'], device_loc, key='DEVICE_ID')
                    device_ids.append(device_id)
                    temps.append(temp)
                logger.debug(f'Processed {len(device_ids)} devices '
                             f'for PC0{series.pcid}')
                time = series['time_recorded'].strftime('%Y-%m-%dT%H:%M:%S%z')
                data.loc[device_ids, 'time_recorded'] = time[:-2]+':'+time[-2:]
                data.loc[device_ids, 'temp'] = temps
                data.loc[device_ids, 'posfid_state'] = series['posfid_state']
    # update temperature value for each device_type for FP temp display
    for mask in [posmask, fidmask]:  # color = temp - temp_mean of device type
        data.loc[mask, 'temp_color'] = (data[mask]['temp']
                                        - data[mask]['temp'].mean(skipna=True))


def process_pc_telemetry_status(query):
    for petal_loc in petal_locs:
        query_petal = query[query['pcid'] == petal_loc]
        if query_petal.empty:
            logger.info(
                f'No status telemetry found for PC0{petal_loc}, resetting...')
            for field in status_colors.keys():
                data_status.loc[field, f'val_{petal_loc}'] = np.nan
                data_status.loc[field, f'color_{petal_loc}'] = 'grey'
        else:
            series = query_petal.iloc[-1]  # select the last row, latest entry
            for field in status_colors.keys():
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
                    logger.error(f'Bad status value: {field} = {value}, '
                                 f'type {type(value)}')
                    continue  # skip bad vaue for this field
                data_status.loc[field, f'val_{petal_loc}'] = series[field]
                data_status.loc[field, f'color_{petal_loc}'] = (
                        status_colors[field][value])


def update_data_and_sources(minutes=5, update_hist=True):
    logger.info(f'Updating plotting sources with new telemetry data '
                f'in the past {minutes} minutes...')
    # do one SQL query for all telemetry data in the past x min for efficiency
    query = query_db(minutes=minutes, table='pc_telemetry_can_all')
    process_pc_telemetry_can_all(query)
    query = query_db(minutes=minutes, table='pc_telemetry_status')
    process_pc_telemetry_status(query)
    source.data = data  # source is a global variable that applies to all plots
    source_status.data = data_status
    if update_hist:
        source_hist.data = generate_hist_data()


def update_plots():
    logger.info(f'Refreshing plots...')
    update_data_and_sources()  # update data using default time span
    logger.debug('Updating plot texts...')
    fp_temp.title.text = (  # update text in plots
        f'Focal Plane Temperature (last updated: {pc.timestamp_str()}, '
        f'refresh interval: {refresh_interval} s)')
    for i, petal_hist in enumerate(petal_hists):
        petal_hist.title.text = (
            f'PC{i:02}, PTL{ptlids[i]} (last updated: {pc.timestamp_str()})')
    logger.info(f'Refreshing in {refresh_interval} s...')


def generate_hist_data(n_bins=10):
    '''create source for histogram, fields are
    left_0_pos, left_0_fid, ..., right_0_pos, ..., top_0_pos,  bottom_0_pos,...
    '''
    data_hist = {}
    for petal_loc in petal_locs:
        petal_data = data[data['petal_loc'] == petal_loc]
        posmask = petal_data['device_type'] == 'POS'
        for device_type, mask in zip(['pos', 'fid'], [posmask, ~posmask]):
            series = petal_data[mask]['temp']
            hist, edges = np.histogram(series[series.notnull()],
                                       density=False, bins=n_bins)
            if np.all(hist == 0):  # empty telemetry data, hist is all zeros
                # use alternative top so bokeh histogram doesn't crash
                hist = hist.astype(np.float64) + 0.1
            data_hist[f'top_{petal_loc}_{device_type}'] = hist
            data_hist[f'left_{petal_loc}_{device_type}'] = edges[:-1]
            data_hist[f'right_{petal_loc}_{device_type}'] = edges[1:]
    return data_hist


def plot_fp_temp():
    tooltips = ([('cursor obsXY', '($x, $y)')]
                + [(col_name, '@'+col_name) for col_name in data.columns
                   if col_name not in ['line_color']])
    fp_temp = figure(title='Focal Plane Temperature',
                     tools='pan,wheel_zoom,reset,hover,save',
                     tooltips=tooltips,
                     aspect_scale=1, plot_width=950, plot_height=1000)
    fp_temp.xaxis.axis_label = 'obsX / mm'
    fp_temp.yaxis.axis_label = 'obsY / mm'
    fp_temp.hover.show_arrow = True
    # low = data['temp_color'].min(skipna=True)
    # high = data['temp_color'].max(skipna=True)
    low, high = -4, 10  # colormap isn't auto-updated when new data come in
    color_mapper = LinearColorMapper(palette=Magma256, low=low, high=high)
    fp_temp.circle(
        x='obs_x', y='obs_y', source=source, radius=5,
        fill_color={'field': 'temp_color', 'transform': color_mapper},
        fill_alpha=0.7, line_color='line_color', line_width=1.8,
        hover_line_color='black')
    colorbar = ColorBar(color_mapper=color_mapper,  # border_line_color=None,
                        ticker=AdaptiveTicker(), orientation='horizontal',
                        title='Deviation from mean of device_type / °C',
                        padding=5, location=(300, 0), height=15, width=250)
    fp_temp.add_layout(colorbar, place='above')  # above
    return fp_temp


def plot_histogram(petal_loc):
    bottom = 0.1  # log cannot properly handle bottom = 0, set it above zero
    p = figure(
        title=f'petal_loc = {petal_loc} Temperature Distribution',
        y_axis_type='log', x_range=(7, 40), y_range=(0.1, 130),
        plot_width=500, plot_height=450)
    for device_type, color in zip(['pos', 'fid'], ['royalblue', 'orangered']):
        p.quad(top=f'top_{petal_loc}_{device_type}', bottom=bottom,
               left=f'left_{petal_loc}_{device_type}',
               right=f'right_{petal_loc}_{device_type}',
               source=source_hist, legend=device_type,
               fill_color=color, line_color="white", alpha=0.5)
    p.y_range.start = bottom
    p.legend.location = "bottom_right"
    p.xaxis.axis_label = 'Temp / °C'
    p.yaxis.axis_label = 'Device Count'
    # add indicator labels only once
    labels = LabelSet(x='x', y='y', text='status_field', source=source_status,
                      x_offset=10, render_mode='css', text_baseline='middle',
                      text_font_size='11pt')
    p.circle(x='x', y='y', radius=0.4, alpha=1, line_alpha=0,
             fill_color=f'color_{petal_loc}', source=source_status)
    p.add_layout(labels)
    return p


# %% main
# note that all variables below are global variables to which funcs have access
pi = PositionerIndex()
conn = psycopg2.connect(host="desi-db", port="5442", database="desi_dev",
                        user="desi_reader", password="reader")  # DB connection
# canonical data (dataframe) and source
logger.info('Initialising data structures...')
data, source, data_status, source_status, posmask, fidmask = init_data()
# load past 5 min and use the latest readout
update_data_and_sources(update_hist=False)
source_hist = ColumnDataSource(generate_hist_data())  # requires updated data
logger.info(f'Initialising plots...')
fp_temp = plot_fp_temp()  # initial plot, focal plane temperature heatmap
# make initial temperature histograms for each petal
petal_hists = [plot_histogram(petal_loc) for petal_loc in petal_locs]
grid = gridplot(petal_hists, ncols=len(petal_locs)//2)
layout = row([fp_temp, grid])
# update plots and save output html and setup curdoc callback
update_plots()
logger.info('Setting up Bokeh document...')
output_file('main.html')
save(layout)
curdoc().add_root(layout)
curdoc().title = 'Focal Plane Telemetry Monitor Application'
# add callback to update existing plots, each webpage creates a callback though
curdoc().add_periodic_callback(update_plots, 1000*refresh_interval)  # in ms
