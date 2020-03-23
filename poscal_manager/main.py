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
from bokeh.plotting import figure
from bokeh.palettes import Magma256, Category10
from bokeh.models import (
    LinearColorMapper, ColorBar, AdaptiveTicker, ColumnDataSource,
    Title, Button, CheckboxButtonGroup)
from bokeh.models.widgets.markups import Div
from bokeh.models.widgets.tables import (
    DataTable, TableColumn, SelectEditor, IntEditor)
from bokeh.layouts import column, layout
import numpy as np
np.seterr(all='raise')


fpa_experts = ['Duan', 'Fanning', 'Fagrelius', 'Schubnell', 'Silber',
               'Poppett', 'Kai']
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
ptls_bt_group = CheckboxButtonGroup(labels=[f'PC{i:02}' for i in range(10)],
                                    active=pcm.pcids, width=300)


def change_selected_calibration(attr, old, new):
    pcm.i_selected = new[0]
    print('selection changed', attr, old, pcm.i_selected)


def on_change_source_data(attr, old, new):
    # old, new, source.data are all the same
    print('Source changed, updating manager data and saving to disk')
    pcm.table = new
    pcm.table_dict_to_df()


def change_ptls(attr, old, new):
    pcm.pcids = new
    print('Changing checked petals', old, pcm.pcids)
    process_calibdf()
    update_histograms_and_scatters()


def plot_heatmap(col):
    '''col is a column name in calibdf, which can be
    R1R2_sum, residuals, gear_ratio_T, gear_ratio_P'''
    # preset data for each quantity to be plotted
    names = {'R1R2_sum': 'R1+R2', 'residuals': 'RMS residuals',
             'GEAR_CALIB_T': 'Gear ratio θ', 'GEAR_CALIB_P': 'Gear ratio φ'}
    units = {'R1R2_sum': ' / mm', 'residuals': ' / mm',
             'GEAR_CALIB_T': '', 'GEAR_CALIB_P': ''}
    lims = {'R1R2_sum': (5.5, 6.5), 'residuals': (0, 0.02),
            'GEAR_CALIB_T': (0.8, 1.2), 'GEAR_CALIB_P': (0.8, 1.2)}
    # begin plot
    data, calibdf = pcm.data, pcm.calibdf
    name, unit, lim = names[col], units[col], lims[col]
    tooltips = ([('cursor obsXY', '($x, $y)')]
                + [(col, '@'+col) for col in calibdf.columns
                   if 'residuals_' not in col])
    heatmap = figure(
        title=f'{name}, expid {data.expid}, {data.mode}',
        tools='pan,wheel_zoom,reset,hover,save', tooltips=tooltips,
        aspect_scale=1, plot_width=450, plot_height=500,
        x_range=(-420, 420), y_range=(-420, 420))
    heatmap.xaxis.axis_label = 'obsX / mm'
    heatmap.yaxis.axis_label = 'obsY / mm'
    heatmap.hover.show_arrow = True
    # low = calibdf[quantity].min(skipna=True)
    # high = calibdf[quantity].max(skipna=True)
    color_mapper = LinearColorMapper(
        palette=Magma256, low=lim[0], high=lim[1])
    heatmap_src = ColumnDataSource(calibdf)
    heatmap.circle(
        x='obs_x', y='obs_y', source=heatmap_src, radius=5,
        fill_color={'field': col, 'transform': color_mapper},
        fill_alpha=0.7, line_color='white', line_width=1.8,
        hover_line_color='black')
    colorbar = ColorBar(
        title=name+unit, color_mapper=color_mapper, ticker=AdaptiveTicker(),
        orientation='horizontal',
        padding=5, location=(0, 0), height=10, width=350)
    heatmap.add_layout(colorbar, place='above')  # above
    return heatmap, heatmap_src


def plot_histogram(col):
    cols = {'R1R2_sum': ['R1R2_sum', 'LENGTH_R1', 'LENGTH_R2'],
            'residuals':  ['residuals', 'residuals_T', 'residuals_P'],
            'GEAR_CALIB_T': ['GEAR_CALIB_T'],
            'GEAR_CALIB_P': ['GEAR_CALIB_P']}
    names = {'R1R2_sum': 'R1+R2', 'residuals': 'RMS residuals',
             'GEAR_CALIB_T': 'Gear ratio θ', 'GEAR_CALIB_P': 'Gear ratio φ'}
    units = {'R1R2_sum': 'mm', 'residuals': 'mm',
             'GEAR_CALIB_T': 'dimensionless', 'GEAR_CALIB_P': 'dimensionless'}
    # lims = {'R1R2_sum': (2.4, 7.2), 'residuals': (0, 50),
    #         'GEAR_CALIB_T': (0.1, 1.4), 'GEAR_CALIB_P': (0.1, 1.4)}
    locs = {'R1R2_sum': (120, 140), 'residuals': 'top_right',
            'GEAR_CALIB_T': 'top_left', 'GEAR_CALIB_P': 'top_left'}
    data, data_hist = pcm.data, pcm.data_hist
    name, unit, loc = names[col], units[col], locs[col]
    tooltips = ([('cursor obsXY', '($x, $y)')]
                + [(k, '@'+k) for k in data_hist.keys()])
    hist = figure(
        title=(f'{name} distribution, expid {data.expid}, {data.mode}, '
               f'{len(pcm.pcids)} petals'),
        tools='pan,wheel_zoom,reset,hover,save', tooltips=tooltips,
        y_axis_type='log', plot_width=450, plot_height=300)  # , x_range=lim)
    hist_src = ColumnDataSource(data_hist)
    colors = iter(Category10[10])
    for i, c in enumerate(cols[col]):
        hist.step(source=hist_src, mode='center', alpha=0.7, legend=c,
                  x=f'x_{c}', y=f'y_{c}',
                  line_color=next(colors))
    hist.xaxis.axis_label = unit
    hist.yaxis.axis_label = 'count / bin'
    hist.legend.location = loc
    hist.legend.padding = 2
    hist.legend.spacing = 1
    return hist, hist_src


def plot_scatter(col):
    names = {'R1R2_sum': 'R1+R2', 'residuals': 'RMS residuals',
             'GEAR_CALIB_T': 'Gear ratio θ', 'GEAR_CALIB_P': 'Gear ratio φ'}
    units = {'R1R2_sum': ' / mm', 'residuals': ' / mm',
             'GEAR_CALIB_T': '', 'GEAR_CALIB_P': ''}
    data, calibdf = pcm.data, pcm.calibdf
    name, unit = names[col], units[col]
    tooltips = ([('cursor obsXY', '($x, $y)')]
                + [(col, '@'+col) for col in calibdf.columns
                   if 'residuals_' not in col])
    scat = figure(
        title=(f'{name} distribution, expid {data.expid}, {data.mode}, '
               f'{len(pcm.pcids)} petals'),
        tools='pan,wheel_zoom,reset,hover,save', tooltips=tooltips,
        plot_width=450, plot_height=300)
    scat_src = ColumnDataSource(calibdf)
    colors = Category10[10]
    for pcid in range(10):  # make 10 lenged items when initialising plot
        scat.circle(x='obs_r', y=f'{col}_{pcid}', source=scat_src, radius=1.5,
                    line_width=0, hover_fill_color='black',
                    color=colors[pcid], legend=f'PC{pcid:02}')
    scat.xaxis.axis_label = 'r / mm'
    scat.yaxis.axis_label = col + unit
    scat.legend.padding = 1
    scat.legend.spacing = 1
    return scat, scat_src


def process_calibdf():
    calibdf = pcm.data.calibdf['FIT'].join(pi)
    calibdf['R1R2_sum'] = calibdf['LENGTH_R1'] + calibdf['LENGTH_R2']
    if 'residuals' not in calibdf.columns:  # arc calibration
        npts = calibdf['residuals_T'][0].size + calibdf['residuals_P'][0].size
        calibdf['residuals'] = np.sqrt(
            ((calibdf['residuals_T']**2).apply(np.nansum)
             + (calibdf['residuals_P']**2).apply(np.nansum))/npts)
    else:  # grid calibration
        calibdf['residuals_T'] = calibdf['residuals_P'] = np.nan
    # for col in ['residuals', 'residuals_T', 'residuals_P']:
    #     calibdf[col] *= 1000  # keep in mm as often times residuals are large
    # more fields for histograms
    n_bins = 40
    data_hist = {}
    mask = calibdf['petal_loc'].isin(pcm.pcids)
    cols = ['R1R2_sum', 'LENGTH_R1', 'LENGTH_R2', 'residuals', 'residuals_T',
            'residuals_P', 'GEAR_CALIB_T', 'GEAR_CALIB_P']
    for col in set(cols) & set(calibdf.columns):
        series = calibdf[mask][col]
        hist, edges = np.histogram(
            series[series.notnull()].apply(np.nanrms),
            density=False, bins=n_bins)
        hist = hist.astype(float)
        hist[hist == 0] += 0.1
        # if np.all(hist == 0):  # hist is all zeros, use alternative top
        #     hist = hist.astype(np.float64) + 0.1  # so quad doesn't crash
        data_hist[f'x_{col}'] = (edges[:-1] + edges[1:])/2
        data_hist[f'y_{col}'] = hist
    # for scatter plots
    calibdf['obs_r'] = np.linalg.norm(calibdf[['obs_x', 'obs_y']], axis=1)
    for col in ({'R1R2_sum', 'residuals', 'GEAR_CALIB_T', 'GEAR_CALIB_P'}
                & set(calibdf.columns)):
        for pcid in range(10):
            mask = calibdf['petal_loc'] == pcid
            if pcid in pcm.pcids:
                calibdf.loc[mask, f'{col}_{pcid}'] = calibdf.loc[mask, col]
            else:
                calibdf.loc[mask, f'{col}_{pcid}'] = np.nan
    pcm.calibdf = calibdf
    pcm.data_hist = data_hist


def update_heatmaps():
    hm_r1r2_src.data = pcm.calibdf
    hm_r1r2.title.text = f'R1+R2, expid {pcm.data.expid}, {pcm.data.mode}'
    hm_res_src.data = pcm.calibdf
    hm_res.title.text = (
        f'RMS residuals, expid {pcm.data.expid}, {pcm.data.mode}')
    if 'arc' in pcm.data.mode:
        print('Updating gear ratio heatmaps for arc calibration')
        hm_grt_src.data = pcm.calibdf
        hm_grt.title.text = (
            f'Gear ratio θ, expid {pcm.data.expid}, {pcm.data.mode}')
        hm_grp_src.data = pcm.calibdf
        hm_grp.title.text = (
            f'Gear ratio φ, expid {pcm.data.expid}, {pcm.data.mode}')


def update_histograms_and_scatters():
    hist_r1r2_src.data = pcm.data_hist
    hist_r1r2.title.text = (
        f'R1+R2 distribution, expid {pcm.data.expid}, {pcm.data.mode}, '
        f'{len(pcm.pcids)} petals')
    scat_r1r2_src.data = pcm.calibdf
    scat_r1r2.title.text = (
        f'R1+R2, expid {pcm.data.expid}, {pcm.data.mode}, '
        f'{len(pcm.pcids)} petals')
    hist_res_src.data = pcm.data_hist
    hist_res.title.text = (
        f'RMS residuals distribution, expid {pcm.data.expid}, '
        f'{pcm.data.mode}, {len(pcm.pcids)} petals')
    scat_res_src.data = pcm.calibdf
    scat_res.title.text = (
        f'RMS residuals, expid {pcm.data.expid}, {pcm.data.mode}, '
        f'{len(pcm.pcids)} petals')
    if 'arc' in pcm.data.mode:
        print('Updating gear ratio histograms for arc calibration')
        hist_grt_src.data = pcm.data_hist
        hist_grt.title.text = (
            f'Gear ratio θ distribution, expid {pcm.data.expid}, '
            f'{pcm.data.mode}, {len(pcm.pcids)} petals')
        scat_grt_src.data = pcm.calibdf
        scat_grt.title.text = (
            f'Gear ratio θ, expid {pcm.data.expid}, {pcm.data.mode}, '
            f'{len(pcm.pcids)} petals')
        hist_grp_src.data = pcm.data_hist
        hist_grp.title.text = (
            f'Gear ratio φ distribution, expid {pcm.data.expid}, '
            f'{pcm.data.mode}, {len(pcm.pcids)} petals')
        scat_grp_src.data = pcm.calibdf
        scat_grp.title.text = (
            f'Gear ratio φ, expid {pcm.data.expid}, {pcm.data.mode}, '
            f'{len(pcm.pcids)} petals')


def update_plots():
    print('current selected', pcm.i_selected, ', pcids', pcm.pcids)
    path = pcm.df.loc[pcm.i_selected, 'data path']
    pcm.data = pcm.read_pickle(path)
    process_calibdf()
    print('Updating plots for', pcm.i_selected, path, pcm.pcids)
    update_heatmaps()
    update_histograms_and_scatters()


path = pcm.df.loc[pcm.i_selected, 'data path']
pcm.data = pcm.read_pickle(path)
process_calibdf()
print('Plotting calibration for', pcm.i_selected, path, pcm.pcids)
# initialise plots, first column, R1+R2
hm_r1r2, hm_r1r2_src = plot_heatmap('R1R2_sum')
hist_r1r2, hist_r1r2_src = plot_histogram('R1R2_sum')
scat_r1r2, scat_r1r2_src = plot_scatter('R1R2_sum')
# second column, fitting residuals, seperated for arc cal and combined for grid
hm_res, hm_res_src = plot_heatmap('residuals')
hist_res, hist_res_src = plot_histogram('residuals')
scat_res, scat_res_src = plot_scatter('residuals')
# third column, gear ratio theta
hm_grt, hm_grt_src = plot_heatmap('GEAR_CALIB_T')
hist_grt, hist_grt_src = plot_histogram('GEAR_CALIB_T')
scat_grt, scat_grt_src = plot_scatter('GEAR_CALIB_T')
# fouth column, gear ratio phi
hm_grp, hm_grp_src = plot_heatmap('GEAR_CALIB_P')
hist_grp, hist_grp_src = plot_histogram('GEAR_CALIB_P')
scat_grp, scat_grp_src = plot_scatter('GEAR_CALIB_P')
# construct webpage layout
col_r1r2 = column([hm_r1r2, hist_r1r2, scat_r1r2])
col_res = column([hm_res, hist_res, scat_res])
col_grt = column([hm_grt, hist_grt, scat_grt])
col_grp = column([hm_grp, hist_grp, scat_grp])
layout = layout([[title],
                 [table],
                 [plot_bt, ptls_bt_group],
                 [col_r1r2, col_res, col_grt, col_grp]])
# add callbacks and layout
source.on_change('data', on_change_source_data)
source.selected.on_change('indices', change_selected_calibration)
plot_bt.on_click(update_plots)
ptls_bt_group.on_change('active', change_ptls)
curdoc().title = 'DESI Positioner Calibration Manager'
curdoc().add_root(layout)
