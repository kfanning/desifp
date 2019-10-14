# -*- coding: utf-8 -*-
'''
Created on Sat Oct 5 09:10:19 2019

@author: Duan Yutong (dyt@physics.bu.edu)

from positioner index and transforms create a df of
    petal_loc, device_loc, device_type, x_dark_1, ..., y_dark_1, ...
device_id

all backlit spots should be close to nominal obsXY.
with every move the spots are still well separated and close to nominal centre.
this doesn't tell if theta/phi motors are sticky or not--it can be done but
would just mean a lot more work. the same thing can be told about the phi motor
when doing 1p calib verification of the 1p movement.
'''


import os
import pandas as pd
from DOSlib.positioner_index import PositionerIndex
from petaltransforms import PetalTransforms

petal_locs = [3]
fvc_img_dir = r'K:\Google Drive\DESI\focal_plane_software\calibration'
bright_exps = ['fvc.20191003213313',  # initial position, close to home
               'fvc.20191003214640',  # After +30 deg phi move
               'fvc.20191003222441']  # After -60 deg theta move
fp_alignment = pd.read_csv(os.path.join(fvc_img_dir,
                                        'focal_plane_alignment_pm.csv'),
                           index_col='PETAL_LOC')
fp_alignment.columns = meapos.columns.str.lower()


def initialise_data():
    path = os.getenv('DOS_POSITIONERINDEXTABLE',
                     '/software/products/PositionerIndexTable-trunk/index_files/desi_positioner_indexes_20190919.csv')
    pi_df = pd.read_csv(path)
    pi_df.columns = pi_df.columns.str.lower()
    # pi_df.set_index('device_id', inplace=True)
    ptlXYZ_df = pd.read_csv(pc.dirs['positioner_locations_file'],
                            usecols=['device_loc', 'X', 'Y', 'Z'],
                            index_col='device_loc')
    data_dfs = []
    for petal_loc in petal_locs:  # create nominal cenre for each petal
        petal_df = (pi_df[pi_df['petal_loc']==petal_loc]
                    .set_index('device_loc').sort_index())
        params = fp_alignment.loc[petal_loc]
        trans = PetalTransforms(Tx=params['XPETAL(mm)'],
                                Ty=params['YPETAL(mm)'],
                                gamma=params['rot_proper(rad)'])
        obsXY = trans.ptlXYZ_to_obsXY(ptlXYZ_df.loc[petal_df.index].T.values)
        petal_df['obsX_0'] = obsXY[0, :]
        petal_df['obsY_0'] = obsXY[1, :]
        data_dfs.append(petal_df)
    data = pd.concat(data_dfs)


def bright_image_spotmatch(bright_exp, i):
    path = os.path.join(fvc_img_dir, bright_exp+'.pos')
    pos = pd.read_csv(path, , index_col=0)  # debug this
    columns = ['x', 'y', 'mag', 'id', 'fwhm']
    pos.set_index('id', inplace=True)
    pxXY = pos[['x', 'y']]
    # for now, assume each fibre spot has a constant id in all exposures
    # need to match id to d
    QS = get_centroids(pxXY)  # centroids in FVC pixel space to FP CS5
    pos['q'], pos['s'] = QS[0, :], QS[1, :]
    obsXY = PetalTransforms.QS_to_obsXY(QS)  # measured 5k obsXY of spots
    pos['obsX'], pos['obsY'] = obsXY[0, :], obsXY[1, :]
    for j, row in pos.iterrows():  # j is id in pos file
        dr = np.linalg.norm(data[['obsX_0', 'obsY_0']].values
                            - row['obsX', 'obsY'].values, axis=1)
        idx = data.iloc[np.argmin(dr)].name  # data index is unique device_id
        data.loc[idx, f'id_{i}'] = j  # spotmatch id
        data.loc[idx, f'obsX_{i}'] = row['obsX']
        data.loc[idx, f'obsY_{i}'] = row['obsY']


def data_consistent():
    '''take average of three pos files, compare with nominal'''
    data['id_consistent'] = (data['id_1'] == data['id_2'] == data['id_3'])
    if np.all(data['id_consistent']):
        print('Spot ids consistent across pos files from bright exposures.')
        return True
    else:
        string = (data[~data['id_consistent']]
                  [['device_id', 'device_loc', 'obsX_0', 'obsY_0',
                    'id_1', 'id_2', 'id_3']])
        print(f'Spot ids inconsistent for:\n{string}')
        return False


def check_movement():
    data['dr_12'] = np.linalg.norm(data[['obsX_1', 'obsY_1']].values
                                   - row['obsX_2', 'obsY_2'].values, axis=1)
    data['dr_23'] = np.linalg.norm(data[['obsX_2', 'obsY_2']].values
                                   - row['obsX_3', 'obsY_3'].values, axis=1)
    data[data['dr_12'] < ]


# %% main
df = initialise_data()
for i, bright_exp in enumrate(bright_exps):
    bright_image_spotmatch(bright_exp, i+1)  # i is the bright exposure index
if data_consistent():
    check_movement

