'''
start server with the following command:
    bokeh serve fp_plots --allow-websocket-origin=desi-2.kpno.noao.edu:5006

view at:
    http://desi-2.kpno.noao.edu:5006/main
'''

import posconstants as pc
import os
from fp_monitor import FPMonitor
from bokeh.io import curdoc  # , output_file, save
from bokeh.layouts import row, gridplot
from bokeh.palettes import Magma256
from bokeh.plotting import figure
from bokeh.models import (
    LinearColorMapper, ColorBar, AdaptiveTicker, LabelSet)


def plot_fp_temp(data, source):
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
    low, high = 15, 30  # colormap isn't auto-updated when new data come in
    color_mapper = LinearColorMapper(palette=Magma256, low=low, high=high)
    fp_temp.circle(
        x='obs_x', y='obs_y', source=source, radius=5,
        fill_color={'field': 'temp_color', 'transform': color_mapper},
        fill_alpha=0.7, line_color='line_color', line_width=1.8,
        hover_line_color='black')
    colorbar = ColorBar(color_mapper=color_mapper,  # border_line_color=None,
                        ticker=AdaptiveTicker(), orientation='horizontal',
                        title='absolute device temperature / °C',
                        padding=5, location=(300, 0), height=15, width=250)
    fp_temp.add_layout(colorbar, place='above')  # above
    return fp_temp


def plot_histogram(source_hist, source_status, petal_loc):
    bottom = 0.1  # log cannot properly handle bottom = 0, set it above zero
    p = figure(
        title=f'petal_loc = {petal_loc} Temperature Distribution',
        y_axis_type='log', x_range=(7, 40), y_range=(0.1, 130),
        plot_width=500, plot_height=450)
    for device_type, color in zip(['pos', 'fid'], ['royalblue', 'orangered']):
        p.quad(top=f'top_{petal_loc}_{device_type}', bottom=bottom,
               left=f'left_{petal_loc}_{device_type}',
               right=f'right_{petal_loc}_{device_type}',
               source=source_hist, legend_label=device_type,
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


def update_plots():
    #print(f'Refreshing plots...')
    if fpm.update_data_and_sources():  # update data using default time span
        #print('Updating plot texts...')
        t_str = pc.timestamp_str(pc.now())
        fp_temp.title.text = (  # update text in plots
            f'Focal Plane Temperature (last updated: {t_str}, '
            f'refresh interval: {refresh_interval} s)')
        for i, petal_hist in enumerate(petal_hists):
            petal_hist.title.text = (
                f'PC{i:02}, PTL{FPMonitor.ptlids[i]} (last updated: {t_str})')
    #print(f'Refreshing in {refresh_interval} s...')


# %% main
refresh_interval = 10
if 'DOS_FPPLOTS_LOGS' not in os.environ:
    os.environ['DOS_FPPLOTS_LOGS'] = os.path.abspath('/data/msdos/desifp/fp_plots')
fpm = FPMonitor(refresh_interval=refresh_interval)
# fpm.logger.info(f'Initialising plots...')
fp_temp = plot_fp_temp(fpm.data, fpm.source)  # initial plot, fp heatmap
# make initial temperature histograms for each petal
petal_hists = [plot_histogram(fpm.source_hist, fpm.source_status, petal_loc)
               for petal_loc in fpm.petal_locs]
grid = gridplot(petal_hists, ncols=len(fpm.petal_locs)//2)
layout = row([fp_temp, grid])
print(f'Updating plots. Refresh interval set to {refresh_interval} s') 
update_plots()
# output_file('main.html')
# save(layout)
curdoc().add_root(layout)
curdoc().title = 'DESI Focal Plane Telemetry Monitor'
# add callback to update existing plots, each webpage creates a callback
curdoc().add_periodic_callback(update_plots, 1000*refresh_interval)
