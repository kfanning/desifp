# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 14:00:33 2019

@author: Duan Yutong
"""


import numpy as np
import matplotlib.pyplot as plt
import posconstants as pc
import petaltransforms
import pandas as pd
trans = petaltransforms.PetalTransforms()
data_ptl = pd.read_csv(pc.dirs['positioner_locations_file'])
data_ptl = data_ptl[['device_loc', 'device_type', 'X', 'Y', 'Z']]
data_flat = data_ptl.copy()
data_flat[['X', 'Y']] = trans.ptlXYZ_to_flatXY(data_ptl[['X', 'Y', 'Z']].values.T).T


def plot_holes(ax, data, color):
    x, y, device_loc, device_type = (data['X'], data['Y'], data['device_loc'],
                                     data['device_type'])
    ax.plot(x, y, '+', color=color, label=device_type.iloc[0])
    for x0, y0, device_loc0, device_type0 in zip(
            x, y, device_loc, device_type):
        circle = plt.Circle((x0, y0), 6.0, ls='--', lw=0.2, color=color,
                            fill=False, label='')
        ax.add_patch(circle)
        ax.annotate(f'{device_loc0}, {device_type0}', xy=(x0-3, y0-3),
                    color=color, size=3)


def plot_map(data, keepout_ptl, keepout_gfa, filename):
    fig, ax = plt.subplots(figsize=(20, 15))
    ax.plot(keepout_ptl[0], keepout_ptl[1], '--', lw=0.3, color='C0',
            marker=',', ms=10, label='KEEPOUT_PTL')
    for x, y in zip(keepout_ptl[0], keepout_ptl[1]):
        ax.annotate(f'({x:.3f}, {y:.3f})', xy=(x, y), color='C0')
    ax.plot(keepout_gfa[0], keepout_gfa[1], '--', lw=0.3, color='C1',
            marker=',', ms=10, label='KEEPOUT_GFA')
    pos = data[data['device_type'] == 'POS']
    fif = data[data['device_type'] == 'FIF']
    gif = data[data['device_type'] == 'GIF']
    opt = data[data['device_type'] == 'OPT']
    etc = data[data['device_type'] == 'ETC']
    for data, color in zip(
            [pos, fif, gif, opt, etc], ['C2', 'C3', 'C4', 'C5', 'C6']):
        plot_holes(ax, data, color)
    # petal below
    ax.plot([20.26, 420], [-0.6, -0.6], '-', lw=0.3, color='C3',
            label='next petal boundary')
    ax.set_aspect('equal')
    ax.legend(loc='center left', prop={'size': 18})
    fig.savefig(filename, bbox_inches='tight')

def add_ptlZ(ptlXY):
    R = np.sqrt(np.square(ptlXY[0, :]) + np.square(ptlXY[1, :]))  # get Z
    Z = pc.R2Z_lookup(R)
    return np.vstack([ptlXY, Z])

if __name__ == '__main__':
    # focalplane/fp_settings/collision_settings/_collision_settings_DEFAULT.conf
    # petal map in ptlXY / obsXY (location 3) coordinates
    KEEPOUT_PTL = np.array(
        [[ 20.260, 410.189, 418.00,  406.3, 399.0, 384.0, 325.837, 20.260,  20.26, 420.0, 420.0, 20.26],
         [  0.000,   0.000,  32.25 ,  89.0, 125.0, 167.5, 235.993, 13.978, 250.00 ,250.0 , -5.0, -5.00]])
    KEEPOUT_GFA = np.array(
        [[295.569, 301.644, 303.538, 305.444, 307.547, 309.204, 320.770, 353.527],
         [207.451, 201.634, 201.588, 201.296, 201.131, 199.968, 184.033, 207.831]])
    plot_map(data_ptl, KEEPOUT_PTL, KEEPOUT_GFA, 'petal_map_ptlxy-obsxy.pdf')
    # petal map in flatXY
    ptlXYZ = add_ptlZ(KEEPOUT_PTL)
    KEEPOUT_PTL_FLAT = trans.ptlXYZ_to_flatXY(ptlXYZ)
    # KEEPOUT_PTL_FLAT[:, 9] = [420, 250]
    ptlXYZ = add_ptlZ(KEEPOUT_GFA)
    KEEPOUT_GFA_FLAT = trans.ptlXYZ_to_flatXY(ptlXYZ)
    plot_map(data_flat, KEEPOUT_PTL_FLAT, KEEPOUT_GFA_FLAT, 'petal_map_flatxy.pdf')
