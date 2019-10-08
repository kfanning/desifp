'''
start server with the following command:
    bokeh serve fp_plots --allow-websocket-origin=desi-2.kpno.noao.edu:5006

view at:
    http://desi-2.kpno.noao.edu:5006/main
'''
from collections import OrderedDict
import os
import numpy as np
import pandas as pd
import posconstants as pc
from petaltransforms import PetalTransforms
import psycopg2
from DOSlib.positioner_index import PositionerIndex
from datetime import datetime, timezone, timedelta
from bokeh.io import curdoc, output_file, save
from bokeh.layouts import row, gridplot
from bokeh.palettes import Magma256
from bokeh.plotting import figure
from bokeh.models import (
    ColumnDataSource, LinearColorMapper, ColorBar, AdaptiveTicker, Label)


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


def strtimenow():
    return datetime.now(timezone.utc).isoformat(timespec=timespec)


def init_data():
    path = os.getenv('DOS_POSITIONERINDEXTABLE',
                     '/software/products/PositionerIndexTable-trunk/index_files/desi_positioner_indexes_20190919.csv')
    pi_df = pd.read_csv(path)
    pi_df.columns = pi_df.columns.str.lower()
    pi_df.set_index('device_id', inplace=True)
    cols = ['spectro_sn', 'spectro_ln',
            'time_recorded', 'temp', 'posfid_state', 'obs_x', 'obs_y',
            'temp_color', 'line_color']
    dtypes = [str, str, 'datetime64[ns]', np.float32, str,
              np.float32, np.float32, np.float32, np.float32, str]
    data = {col: pd.Series(dtype=dt) for col, dt in zip(cols, dtypes)}
    data = pi_df.join(pd.DataFrame(data=data, index=pi_df.index))
    posmask = data['device_type'] == 'POS'
    data.loc[posmask, 'line_color'] = 'white'
    data.loc[~posmask, 'line_color'] = 'green'
    # add obsXYZ positions to data for plotting
    path = os.path.join(os.getenv('PLATE_CONTROL_DIR',
                                  '/software/products/plate_control-trunk'),
                        'petal', 'positioner_locations_0530v14.csv')
    ptlXYZ_df = pd.read_csv(path,
                            usecols=['device_location_id', 'X', 'Y', 'Z'],
                            index_col='device_location_id')
    ptlXYZ_df.index.rename('device_loc', inplace=True)
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
    # initial pc telemetry status data
    data_status = {field: np.nan for field in
                   ['time_recorded'] + list(status_colors.keys())}
    for i, field in enumerate(status_colors.keys()):
        data_status[field+'_color'] = 'grey'  # initial colors
        data_status[field+'_x'] = 23 + 12 * (i//9)
        data_status[field+'_y'] = np.power(10, 2 - 0.17 * (i % 9))
    sources_status = {petal_loc: data_status for petal_loc in petal_locs}
    return data, ColumnDataSource(data), sources_status


def query_db(hours=24, table='pc_telemetry_can_all'):
    time_range = timedelta(hours=hours)
    time_cut = (datetime.utcnow() - time_range).strftime('%Y-%m-%d %H:%M:%S')
    query = pd.read_sql_query(
        f"SELECT * FROM {table} WHERE time_recorded >= '{time_cut}'",
        conn)
    return query_db(hours=hours+1, table=table) if query.empty else query


def process_pc_telemetry_can_all(query):
    '''adds the result of a sql query to the data source'''
    for i, series in query.iterrows():
        device_ids, temps = [], []
        for device_loc, temp in series['posfid_temps'].items():
            if 'can' in device_loc:
                print(f"Skipping DB entry in old telemetry format "
                      f"submitted by PC-{series['pcid']} "
                      f"at {series['time_recorded']}")
                break
            device_ids.append(pi.find_by_petal_loc_device_loc(
                series['pcid'], device_loc, key='DEVICE_ID'))
            temps.append(temp)
        time = series['time_recorded'].strftime('%Y-%m-%dT%H:%M:%S%z')
        data.loc[device_ids, 'time_recorded'] = time[:-2] + ':' + time[-2:]
        data.loc[device_ids, 'temp'] = temps
        data.loc[device_ids, 'posfid_state'] = series['posfid_state']
    # update mean temperature of each device_type
    posmask = data['device_type'] == 'POS'
    for mask in [posmask, ~posmask]:  # color = temp / temp_mean of device type
        data.loc[mask, 'temp_color'] = (data[mask]['temp']
                                        - data[mask]['temp'].mean(skipna=True))


def process_pc_telemetry_status(query):
    for petal_loc in petal_locs:
        query_petal = query[query['pcid'] == petal_loc]
        if query_petal.empty:
            continue  # no data for this petal/petalcontroller, skip to next
        # select the last row, latest entry
        series = query_petal.sort_values('time_recorded').iloc[-1]
        for field in ['time_recorded'] + list(status_colors.keys()):
            sources_status[petal_loc][field] = series[field]
        for field in status_colors.keys():
            value = sources_status[petal_loc][field]
            if not (value == 0 or value == 1 or np.isnan(value)):
                print(f'Bad status value: '
                      f'{field} = {value}, type {type(value)}')
                continue
            sources_status[petal_loc][field+'_color'] = (
                    status_colors[field][value])


def update_data_and_source(hours=24, update_hist=True):
    query = query_db(hours=hours, table='pc_telemetry_can_all')
    process_pc_telemetry_can_all(query)
    query = query_db(hours=hours, table='pc_telemetry_status')
    process_pc_telemetry_status(query)
    source.data = data  # source is a global variable that applies to all plots
    if update_hist:
        source_hist.data = generate_hist_data()


def update_plots():  # minute
    print('Refreshing plots...')
    update_data_and_source(hours=0.1)  # update data
    plot_status_indicators()
    fp_temp.title.text = (  # update text in plots
        f'Focal Plane Temperature (last updated: {strtimenow()}, '
        f'refresh interval: {refresh_interval} s)')
    for i, petal_hist in enumerate(petal_hists):
        petal_hist.title.text = (
            f'petal_loc = {i} (last updated: {strtimenow()})')
    print(f'Refreshing in {refresh_interval} s...')


def generate_hist_data(n_bins=15):
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
            data_hist[f'top_{petal_loc}_{device_type}'] = hist
            data_hist[f'left_{petal_loc}_{device_type}'] = edges[:-1]
            data_hist[f'right_{petal_loc}_{device_type}'] = edges[1:]
    return data_hist


def plot_histogram(petal_loc):
    bottom = 0.1  # log cannot properly handle bottom = 0, set above zero
    p = figure(
        title=f'petal_loc = {petal_loc} Temperature Distribution',
        y_axis_type='log', x_range=(9, 45), y_range=(0.1, 130),
        plot_width=500, plot_height=450)
    for device_type, color in zip(['pos', 'fid'], ['royalblue', 'orangered']):
        p.quad(top=f'top_{petal_loc}_{device_type}',
               bottom=bottom,
               left=f'left_{petal_loc}_{device_type}',
               right=f'right_{petal_loc}_{device_type}',
               source=source_hist, legend=device_type,
               fill_color=color, line_color="white", alpha=0.5)
    p.y_range.start = bottom
    p.legend.location = "top_left"
    p.xaxis.axis_label = 'Temp / °C'
    p.yaxis.axis_label = 'Device Count'
    # add indicator labels only once
    for field in status_colors.keys():
        label = Label(x=sources_status[petal_loc][field+'_x']+1,
                      y=sources_status[petal_loc][field+'_y'],
                      text=field, render_mode='css', text_baseline='middle',
                      text_font_size='11pt')
        p.add_layout(label)
    return p


def plot_status_indicators():
    for petal_loc, p in enumerate(petal_hists):
        for field in status_colors.keys():
            p.circle(x=sources_status[petal_loc][field+'_x'],
                     y=sources_status[petal_loc][field+'_y'],
                     fill_color=sources_status[petal_loc][field+'_color'],
                     line_alpha=0,
                     radius=0.5, alpha=1)


# %% main
# note that all variables below are global variables to which funcs have access
pi = PositionerIndex()
conn = psycopg2.connect(host="desi-db", port="5442", database="desi_dev",
                        user="desi_reader", password="reader")  # DB connection
# canonical data (dataframe) and source
data, source, sources_status = init_data()
# load past 24 hr and use the latest readout
update_data_and_source(hours=24, update_hist=False)
source_hist = ColumnDataSource(generate_hist_data())  # requires updated data
# initial plot, focal plane temperature heatmap
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
color_mapper = LinearColorMapper(palette=Magma256,
                                 low=data['temp_color'].min(skipna=True),
                                 high=data['temp_color'].max(skipna=True))
fp_temp.circle(
    x='obs_x', y='obs_y', source=source, radius=5,
    fill_color={'field': 'temp_color', 'transform': color_mapper},
    fill_alpha=0.7, line_color='line_color', line_width=1.3,
    hover_line_color='black')
colorbar = ColorBar(color_mapper=color_mapper,  # border_line_color=None,
                    ticker=AdaptiveTicker(), orientation='horizontal',
                    title='Deviation from mean of device_type / °C',
                    padding=5, location=(300, 0), height=15, width=250)
fp_temp.add_layout(colorbar, place='above')  # above
# make temperature histograms for each petal
petal_hists = [plot_histogram(petal_loc) for petal_loc in petal_locs]
plot_status_indicators()
grid = gridplot(petal_hists, ncols=len(petal_locs)//2)
layout = row([fp_temp, grid])
update_plots()
output_file('main.html')
save(layout)
curdoc().add_root(layout)
curdoc().title = 'Focal Plane Telemetry Monitor Application'
# add callback to update existing plots
curdoc().add_periodic_callback(update_plots, 1000*refresh_interval)  # in ms
