# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 15:14:18 2020

@author: Duan Yutong (dyt@physics.bu.edu)

start server with the following command:

bokeh serve poscal_manager --allow-websocket-origin=desi-2.kpno.noao.edu:5006

view at: http://desi-2.kpno.noao.edu:5006/poscal_manager
"""

import poscal_manager as pm
from bokeh.io import curdoc  # , output_file, save
# from bokeh.layouts import row
from bokeh.plotting import figure
from bokeh.palettes import Magma256
from bokeh.models import (
    LinearColorMapper, ColorBar, AdaptiveTicker, LabelSet)
from bokeh.models import ColumnDataSource, Button, CheckboxButtonGroup
from bokeh.models.widgets.markups import Div
from bokeh.models.widgets.tables import (
    DataTable, TableColumn, SelectEditor, IntEditor)
from bokeh.layouts import layout
import numpy as np
np.seterr(all='raise')


pcids = list(range(10))
fpa_experts = ['Duan', 'Fanning', 'Fagrelius', 'Schubnell',
               'Silber', 'Poppett']
pcm = pm.PosCalManager()
pi = pm.get_positioner_index()
source = ColumnDataSource(data=pcm.table)
source.selected.indices = [pcm.i_selected]  # ['1d']['indices'] = [i_selected]
title = Div(text='''
<font size="4">Positioner Calibrations</font> (some columns editable)''',
            width=500)
columns = [TableColumn(field='UTC', title='UTC', width=160),
           TableColumn(field='expid', title='expid', width=50),
           TableColumn(field='test name', title='test name', width=260),
           TableColumn(field='dome', title='dome', width=50,
                       editor=SelectEditor(options=['open', 'closed', '?'])),
           TableColumn(field='zenith angle', title='zenith angle', width=69,
                       editor=IntEditor()),
           TableColumn(field='tracking', title='tracking', width=50,
                       editor=SelectEditor(options=['on', 'off', '?'])),
           TableColumn(field='PMTC', title='PMTC', width=50,
                       editor=SelectEditor(options=['on', 'off', '?'])),
           TableColumn(field='PCVF', title='PCVF', width=50,
                       editor=SelectEditor(options=['on', 'off', '?'])),
           TableColumn(field='dome louvers', title='dome louvers', width=75,
                       editor=SelectEditor(options=['on', 'off', '?'])),
           TableColumn(field='dome fans', title='dome fans', width=60,
                       editor=SelectEditor(options=['on', 'off', '?'])),
           TableColumn(field='B29 fan', title='B29 fan', width=50,
                       editor=SelectEditor(options=['on', 'off', '?'])),
           TableColumn(field='cage baffle', title='cage baffle', width=63,
                       editor=SelectEditor(options=['on', 'off', '?'])),
           TableColumn(field='barrel insulation', title='barrel insulation',
                       width=85,
                       editor=SelectEditor(options=['on', 'off', '?'])),
           TableColumn(field='operator', title='operator',
                       width=70,
                       editor=SelectEditor(options=sorted(fpa_experts)+['?'])),
           TableColumn(field='data path', title='data path', width=700)]
table = DataTable(source=source, columns=columns, editable=True,
                  sortable=False, reorderable=False, fit_columns=False,
                  default_size=1300, height=min(27*pcm.len+30, 600),
                  min_width=1300, sizing_mode='stretch_width')
plot_bt = Button(label="Plot Selected Calibration", button_type='primary',
                 width=300)
ptls_bt_group = CheckboxButtonGroup(labels=[f'PC{i:02}' for i in pcids],
                                    active=pcids, width=300)


def change_selected_calibration(attr, old, new):
    pcm.i_selected = new[0]
    print('selection changed', attr, old, pcm.i_selected)


def on_change_source_data(attr, old, new):
    # old, new, source.data are all the same
    print('Source changed, updating manager data and saving to disk')
    pcm.table = new
    pcm.table_dict_to_df()


def change_ptls(attr, old, new):
    pcids = new
    print('Changing checked petals', old, pcids)


names = {'R1R2_sum': 'R1+R2', 'residuals': 'RMS residuals',
         'GEAR_CALIB_T': 'Gear ratio θ', 'GEAR_CALIB_P': 'Gear ratio φ'}
units = {'R1R2_sum': ' / mm', 'residuals': ' / μm',
         'GEAR_CALIB_T': '', 'GEAR_CALIB_P': ''}
lims = {'R1R2_sum': (5.5, 6.5), 'residuals': (0, 20),
        'GEAR_CALIB_T': (0.8, 1.2), 'GEAR_CALIB_P': (0.8, 1.2)}


def plot_heatmap(data, col):
    '''col is a column name in calibdf, which can be
    R1R2_sum, residuals, gear_ratio_T, gear_ratio_P'''
    name, unit, lim = names[col], units[col], lims[col]
    calibdf = process_calibdf(data.calibdf['FIT'])
    calibdf['residuals'] *= 1000
    tooltips = ([('cursor obsXY', '($x, $y)')]
                + [(col, '@'+col) for col in calibdf.columns
                   if 'residuals_' not in col])
    title = f'{name}, expid {data.expid}, {data.mode} calibration'
    heatmap = figure(
        title=title, tools='pan,wheel_zoom,reset,hover,save',
        x_range=(-420, 420), y_range=(-420, 420), tooltips=tooltips,
        aspect_scale=1, plot_width=450, plot_height=500)
    heatmap.xaxis.axis_label = 'obsX / mm'
    heatmap.yaxis.axis_label = 'obsY / mm'
    heatmap.hover.show_arrow = True
    # low = calibdf[quantity].min(skipna=True)
    # high = calibdf[quantity].max(skipna=True)
    color_mapper = LinearColorMapper(
        palette=Magma256, low=lim[0], high=lim[1])
    heatmap_source = ColumnDataSource(calibdf)
    heatmap.circle(
        x='obs_x', y='obs_y', source=heatmap_source, radius=5,
        fill_color={'field': col, 'transform': color_mapper},
        fill_alpha=0.7, line_color='white', line_width=1.8,
        hover_line_color='black')
    colorbar = ColorBar(
        title=name+unit, color_mapper=color_mapper,
        ticker=AdaptiveTicker(), orientation='horizontal',
        padding=5, location=(0, 0), height=10, width=330)
    heatmap.add_layout(colorbar, place='above')  # above
    return heatmap, heatmap_source


def plot_histogram(data, quantity):
    bottom = 0.1  # log cannot properly handle bottom = 0, set it above zero
    hist = figure(
        title=f'{quantity} distribution, {len(pcid)} petals',
        y_axis_type='log', plot_width=450, plot_height=250)
    # fields = [quantity] if 
    # for device_type, color in zip(['pos', 'fid'], ['royalblue', 'orangered']):
    #     p.quad(top=f'top_{petal_loc}_{device_type}', bottom=bottom,
    #            left=f'left_{petal_loc}_{device_type}',
    #            right=f'right_{petal_loc}_{device_type}',
    #            source=source_hist, legend=device_type,
    #            fill_color=color, line_color="white", alpha=0.5)
    # p.y_range.start = bottom
    # p.legend.location = "bottom_right"
    # p.xaxis.axis_label = 'Temp / °C'
    # p.yaxis.axis_label = 'Device Count'
    # # add indicator labels only once
    # labels = LabelSet(x='x', y='y', text='status_field', source=source_status,
    #                   x_offset=10, render_mode='css', text_baseline='middle',
    #                   text_font_size='11pt')
    # p.circle(x='x', y='y', radius=0.4, alpha=1, line_alpha=0,
    #          fill_color=f'color_{petal_loc}', source=source_status)
    # p.add_layout(labels)
    # return p


def process_calibdf(calibdf):
    calibdf = calibdf.join(pi)
    calibdf['R1R2_sum'] = calibdf['LENGTH_R1'] + calibdf['LENGTH_R2']
    if 'residuals' not in calibdf.columns:
        npts = calibdf['residuals_T'][0].size + calibdf['residuals_P'][0].size
        calibdf['residuals'] = np.sqrt(
            ((calibdf['residuals_T']**2).apply(np.nansum)
             + (calibdf['residuals_P']**2).apply(np.nansum))/npts)
    return calibdf


def update_plots():
    print('current selected', pcm.i_selected)
    path = pcm.df.loc[pcm.i_selected, 'data path']
    data = pcm.read_pickle(path)
    calibdf = process_calibdf(data.calibdf['FIT'])
    print('Updating plots for', pcm.i_selected, path, pcids)
    hm_r2r2_src.data = calibdf
    hm_r1r2.title.text = f'R1+R2, expid {data.expid}, {data.mode} calibration'
    hm_res_src.data = calibdf
    hm_res.title.text = (
        f'RMS residuals, expid {data.expid}, {data.mode} calibration')
    if 'arc' in data.mode:
        print('Updating gear ratio plots for arc calibration')
        hm_grt_src.data = calibdf
        hm_grt.title.text = (
            f'Gear ratio θ, expid {data.expid}, {data.mode} calibration')
        hm_grp_src.data = calibdf
        hm_grp.title.text = (
            f'Gear ratio φ, expid {data.expid}, {data.mode} calibration')


path = pcm.df.loc[pcm.i_selected, 'data path']
print('Plotting calibration for', pcm.i_selected, path, pcids)
hm_r1r2, hm_r2r2_src = plot_heatmap(pcm.read_pickle(path), 'R1R2_sum')
hm_res, hm_res_src = plot_heatmap(pcm.read_pickle(path), 'residuals')
hm_grt, hm_grt_src = plot_heatmap(pcm.read_pickle(path), 'GEAR_CALIB_T')
hm_grp, hm_grp_src = plot_heatmap(pcm.read_pickle(path), 'GEAR_CALIB_P')
source.on_change('data', on_change_source_data)
source.selected.on_change('indices', change_selected_calibration)
plot_bt.on_click(update_plots)
ptls_bt_group.on_change('active', change_ptls)
layout = layout([[title],
                 [table],
                 [plot_bt, ptls_bt_group],
                 [hm_r1r2, hm_res, hm_grt, hm_grp]])
curdoc().title = 'DESI Positioner Calibration Manager'
curdoc().add_root(layout)
# output_file('main.html')
# save(table)
