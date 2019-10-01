# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 20:49:23 2019

@author: givoltage
"""
import os
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from petaltransforms import PetalTransforms
# sample fid position from Klaus
alignment_kent = {  # PTL11, location 6
    "petal_id": 11,
    "petal_offset_x": 0.013,
    "petal_offset_y": 0.065,
    "petal_offset_z": 0.0,
    "petal_rot_1": 0.0,
    "petal_rot_2": 0.0,
    "petal_rot_3": 1.8844738812803254}
fid_centre_def = 'pinhole4'  # 'center'
data_dir = r'/home/msdos/fp_temp_files/'
# data_dir = r'Downloads'
data = pd.read_json(os.path.join(data_dir, 'petal11.json'), orient='index')
data['id'] = 6000 + data['device_loc']
data = data.reset_index().set_index('id')
data_old = (pd.read_csv(os.path.join(data_dir, 'fiducial-desi4.dat'),
                        sep=r'\s+')
            .rename(columns={'#id': 'id'}).set_index('id'))
data_kent = (pd.read_csv(os.path.join(data_dir, 'fiducial-desi4_svn.dat'),
                         sep=r'\s+', skiprows=[1])
             .rename(columns={'#SERIAL': 'id'}).set_index('id'))
data['x'] = [d['x'] for d in data[fid_centre_def]]
data['y'] = [d['y'] for d in data[fid_centre_def]]
data['z'] = [d['z'] for d in data[fid_centre_def]]
data['q_old'] = data_old['q']
data['s_old'] = data_old['s']
data['q_kent'] = data_kent['Q']
data['s_kent'] = data_kent['S']


def calculate_qst(alignment):
    trans = PetalTransforms(Tx=alignment[0],
                            Ty=alignment[1],
                            Tz=alignment[2],
                            alpha=alignment[3],
                            beta=alignment[4],
                            gamma=alignment[5])
    obsXYZ = trans.ptlXYZ_to_obsXYZ(data[['x', 'y', 'z']].values.T)
    data['obs_x'] = obsXYZ[0, :]
    data['obs_y'] = obsXYZ[1, :]
    data['obs_z'] = obsXYZ[2, :]
    QST = trans.obsXYZ_to_QST(obsXYZ)
    data['qst_q'] = QST[0, :]
    data['qst_s'] = QST[1, :]
    data['qst_t'] = QST[2, :]
    QS = trans.ptlXYZ_to_QS(data[['x', 'y', 'z']].values.T)
    data['q_new'], data['s_new'] = QS[0, :], QS[1, :]
    data['dq_old_new'] = data['q_new'] - data['q_old']
    data['ds_old_new'] = data['s_new'] - data['s_old']
    data['dq_old'] = data['q_old'] - data['q_kent']
    data['ds_old'] = data['s_old'] - data['s_kent']
    data['dq_new'] = data['q_new'] - data['q_kent']
    data['ds_new'] = data['s_new'] - data['s_kent']
    return np.sqrt(np.mean(np.square(data[['dq_new', 'ds_new']].values)))


# minimise
p0 = np.array([0.013, 0.065, 0.0, 0.0, 0.0, 1.8844738812803254])
calculate_qst(p0)
data.to_csv('qst_test.csv')

bounds = ((-10, 10), (-10, 10), (-10, 10),
          (-np.pi/2, np.pi//2), (-np.pi/2, np.pi/2), (-np.pi, np.pi))
solution = minimize(calculate_qst, p0, bounds=bounds, method='SLSQP',
                    options={'disp': True, 'maxiter': 1000})
print(f'Best-fit transformation params:\n{solution.x}')
calculate_qst(solution.x)
data.to_csv('qst_test_minimised.csv')
