# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 15:40:51 2019

@author: givoltage
"""

import os
import numpy as np
import pandas as pd
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.colors as colors

missing_pos = [['M02956', 'M04529', 'M06783'],
               ['M03156', 'M04133'],
               ['M03120', 'M04283', 'M05607'],
               ['M03294', 'M03757', 'M03966', 'M03990', 'M04243', 'M04478'],
               ['M01190', 'M02792', 'M03493', 'M03596', 'M03599', 'M05152'],
               ['M01727', 'M06214', 'M07163'],
               ['M03767', 'M05638', 'M06009', 'M06847'],
               ['M02582', 'M03159', 'M03554', 'M03837', 'M03943', 'M04217'],
               ['M02732', 'M04583'],
               ['M05070', 'M07601']]
overexposed = ['M03554', 'M05152', 'M07163']
pos_filename = 'fvc.20191004001051.pos'  # this is when the bright image was taken
PC00 = ['fvc.20191003232306.pos', 'fvc.20191003233356.pos', 'fvc.20191003234547.pos']
PC01 = ['fvc.20191003234904.pos', 'fvc.20191004000051.pos', 'fvc.20191004001051.pos']
PC02 = ['fvc.20191004170751.pos', 'fvc.20191004190200.pos', 'fvc.20191004191322.pos']
PC03 = ['fvc.20191003213313.pos', 'fvc.20191003214640.pos', 'fvc.20191003222441.pos']
PC04 = ['fvc.20191003222441.pos', 'fvc.20191004190200.pos', 'fvc.20191007152624.pos']
PC05 = ['fvc.20191007152624.pos', 'fvc.20191007153833.pos', 'fvc.20191007162144.pos']
PC06 = ['fvc.20191007162144.pos', 'fvc.20191007172910.pos', 'fvc.20191007180412.pos']
PC07 = ['fvc.20191007180412.pos', 'fvc.20191007183104.pos', 'fvc.20191007184307.pos']
PC08 = ['fvc.20191007184307.pos', 'fvc.20191007185532.pos', 'fvc.20191007190739.pos']
PC09 = ['fvc.20191007190739.pos', 'fvc.20191007192116.pos', 'fvc.20191007193250.pos']
petlist = [PC00, PC01, PC02, PC03, PC04, PC05, PC06, PC07, PC08, PC09]
petal_locs = range(10)
path = r'K:\Google Drive\DESI\focal_plane_software\bright image calibration\fvc.20191004115421.fits'
# FVC image has a y-flip w.r.t. y-axis compared to CS5
image_data = fits.getdata(path, ext=0)
data = pd.read_pickle('data_spotmatch.pkl', compression='gzip')


# Define a triangle that encloses each petal in FVC pixels
def getpostriangle(petal_loc):
    xc, yc = 3020.3600823, 3015
    theta0 = 180.4
    theta1 = (36*petal_loc - 18. + theta0)*np.pi/180.
    theta2 = (36*petal_loc + 18. + theta0)*np.pi/180.
    x1 = 3000*np.sin(theta1) + xc
    y1 = 3000*np.cos(theta1) + yc
    x2 = 3000*np.sin(theta2) + xc
    y2 = 3000*np.cos(theta2) + yc
    return [x1, x2, xc, x1], [y1, y2, yc, y1]


def determine_exp_no(petal_loc):
    number = int(pos_filename[4:-4])
    numbers = np.array([int(fn[4:-4]) for fn in petlist[petal_loc]])
    if np.all(numbers < number):
        # all exposures were taken before the bright image
        exp_no = 2
    elif np.all(numbers > number):
        # all exposures were taken after the bright image
        exp_no = 0
    else:
        exp_no = int(np.where((numbers <= number) == True)[0][-1])
    return exp_no


def plot_dark_pos(posid):
    ext_margin = 200
    petal_loc = data.loc[posid, 'petal_loc']
    xcentre = data.loc[posid, 'obs_x_px']
    ycentre = data.loc[posid, 'obs_y_px']
    xmin, xmax = xcentre-ext_margin, xcentre+ext_margin
    ymin, ymax = ycentre-ext_margin, ycentre+ext_margin
    if posid in overexposed:
        vmin, vmax = 2e3, 5e3
    else:
        vmin, vmax = 1e3, 2.5e3
    fig, ax = plt.subplots(figsize=(16, 16))
    ax.imshow(image_data, cmap='gray', aspect='equal', origin='lower',
              norm=colors.LogNorm(vmin=vmin, vmax=vmax))
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])
    # fig.colorbar(img)
    for petal_loc in petal_locs:
        xco, yco = getpostriangle(petal_loc)
        ax.plot(xco, yco, 'y--')
        dp = data[data['petal_loc'] == petal_loc]
        exp_no = determine_exp_no(petal_loc)
        ax.plot(xcentre, ycentre, 'ro', fillstyle='none', ms=200, mew=1)
        ax.plot(dp[f'fvc_x_{petal_loc}{exp_no}'],
                dp[f'fvc_y_{petal_loc}{exp_no}'],
                'go', fillstyle='none', ms=15, mew=3)
    return fig, ax


if __name__ == '__main__':
    for petal_loc in petal_locs:
        for posid in missing_pos[petal_loc]:
            if posid not in overexposed:
                continue
            print(f'Click on estimated fibre position for {posid}')
            fig, ax = plot_dark_pos(posid)
            # get one point, the 0th point in the list
            x, y = plt.ginput(1, timeout=0)[0]
            data.loc[posid, 'fvc_x_click'] = x
            data.loc[posid, 'fvc_y_click'] = y
            print(f'FVC pixel position ({x}, {y}) recorded for {posid}')
            ax.plot(x, y, 'rx')
            plt.draw()
            plt.show()  # show updated canvas draw, user to close plot window
            fig.savefig(os.path.join('dark_pos',
                                     f'{posid}-{x:.3f}_{y:.3f}.png'),
                        bbox_inches='tight', dpi=300)
            plt.close('all')
    data.to_csv('data_spotmatch_click.csv')
    data.to_pickle('data_spotmatch_click.pkl', compression='gzip')
