'''
start server with the following command:
    bokeh serve --allow-websocket-origin=desi-2.kpno.noao.edu:5006 main.py

view at:
    http://desi-2.kpno.noao.edu:5006/main
'''
import os
import numpy as np
import pandas as pd
import posconstants as pc
from petaltransforms import PetalTransforms
import psycopg2
from DOSlib.positioner_index import PositionerIndex
from datetime import datetime, timezone, timedelta
# from bokeh.io import show
from bokeh.io import curdoc, show
from bokeh.layouts import row, column
from bokeh.palettes import Magma256
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, LinearColorMapper


def init_data():
    filepath = os.getenv('DOS_POSITIONERINDEXTABLE',
                         '/software/products/PositionerIndexTable-trunk/index_files/desi_positioner_indexes_20190919.csv')
    pi_df = pd.read_csv(filepath)
    pi_df.columns = pi_df.columns.str.lower()
    pi_df.set_index('device_id', inplace=True)
    cols = ['time', 'temp', 'posfid_state', 'obs_x', 'obs_y', 'temp_color']
    dtypes = ['datetime64[ns]', np.float32, str,
              np.float32, np.float32, np.float32, np.float32]
    data = {col: pd.Series(dtype=dt) for col, dt in zip(cols, dtypes)}
    data = pi_df.join(pd.DataFrame(data=data, index=pi_df.index))
    # add obsXYZ positions to data for plotting
    df = pd.read_csv(pc.dirs['positioner_locations_file'],
                     usecols=['device_loc', 'X', 'Y', 'Z'],
                     index_col='device_loc')
    for petal_loc in range(10):
        trans = PetalTransforms(gamma=np.pi/5*(petal_loc-3))
        obsXY = trans.ptlXYZ_to_obsXYZ(df.T.values)[:2, :]
        xy_df = pd.DataFrame(data=obsXY.T, index=df.index,
                             columns=['obs_x', 'obs_y'])
        petal_data = (data[data['petal_loc'] == petal_loc]
                      .reset_index().set_index('device_loc'))
        petal_data.update(xy_df)
        data.update(petal_data.reset_index().set_index('device_id'))
    return data, ColumnDataSource(data)


def query_db(hours=24):
    time_range = timedelta(hours=hours)
    time_cut = (datetime.utcnow()-time_range).strftime("%Y-%m-%d %H:%M:%S")
    query = pd.read_sql_query(
        f"""SELECT * FROM pc_telemetry_can_all
            WHERE time_recorded >= '{time_cut}'""",
        conn)
    if query.empty:
        return query_db(hours=hours+1)
    return query


def process_query_data(query):
    for i, series in query.iterrows():
        device_ids, temps = [], []
        for device_loc, temp in series['posfid_temps'].items():
            if 'can' in device_loc:
                print(f"Skipping entry at {series['time_recorded']}"
                      " in old telemetry format")
                break
            device_ids.append(pi.find_by_petal_loc_device_loc(
                series['pcid'], device_loc, key='DEVICE_ID'))
            temps.append(temp)
        data.loc[device_ids, 'time'] = series['time_recorded'].isoformat()
        data.loc[device_ids, 'temp'] = temps
        data.loc[device_ids, 'posfid_state'] = series['posfid_state']
    # update mean temperature of each device_type
    posmask = data['device_type'] == 'POS'
    for mask in [posmask, ~posmask]:  # color = temp / temp_mean of device type
        data.loc[mask, 'temp_color'] = (data[mask]['temp']
                                        / data[mask]['temp'].mean(skipna=True))


def update_data_and_source(hours=24):
    query = query_db(hours=hours)
    process_query_data(query)
    source.data = data


def update_plots():  # minute
    print('Refreshing plots...')
    update_data_and_source(hours=0.1)
    fp_temp.title.text = (
        'Focal Plane Temperature\n'
        f' (last updated: {datetime.now(timezone.utc).isoformat()})')
    print(f'Refreshing in {refresh_interval} s...')


# %% main
refresh_interval = 10  # seconds
pi = PositionerIndex()
conn = psycopg2.connect(host="desi-db", port="5442", database="desi_dev",
                        user="desi_reader", password="reader")  # DB connection
data, source = init_data()  # data and source contain empty telemetry data now
update_data_and_source(hours=24)  # load past 24 hr and use the latest readout
# initial plot, focal plane temperature heatmap
tools = 'pan,wheel_zoom,reset,hover,save'
title = 'Focal Plane Temperature'
tooltips = ([('obsXY', '($x, $y)')]
            + [(col_name, '@'+col_name) for col_name in data.columns])
fp_temp = figure(title=title, tools=tools, tooltips=tooltips,
                 plot_width=1000, plot_height=1000)
fp_temp.hover.show_arrow = True
# fp_temp.hover.point_policy = "follow_mouse"
color_mapper = LinearColorMapper(palette=Magma256)
fp_temp.circle(
    x='obs_x', y='obs_y', source=source, radius=5,
    fill_color={'field': 'temp_color', 'transform': color_mapper},
    fill_alpha=0.7, line_color='white', line_width=0.5,
    hover_line_color='black')
layout = row(fp_temp)
curdoc().add_root(layout)
curdoc().title = title
# show(fp_temp)
curdoc().add_periodic_callback(update_plots, 1000*refresh_interval)
