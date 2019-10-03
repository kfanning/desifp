# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 20:49:23 2019

@author: Duan Yutong (dyt@physics.bu.edu)
"""
import os
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from petaltransforms import PetalTransforms

fid_centre_def = 'center'
data_dir = r'/home/msdos/fp_temp_files/'
# data_dir = r'Downloads'
petal_locs = {2: 7, 3: 3, 4: 0, 5: 1, 6: 2, 7: 8, 8: 4, 9: 9, 10: 5, 11: 6}
data_kent = (pd.read_csv(os.path.join(data_dir, 'fiducial-desi.dat'),
                         sep=r'\s+', skiprows=[1])
             .rename(columns={'#SERIAL': 'id'}).set_index('id'))

def qst_petal(petal_id):
    data = pd.read_json(os.path.join(data_dir, f'petal{petal_id}.json'),
                        orient='index')
    petal_loc = petal_locs[petal_id]
    data['id'] = petal_loc*1000 + data['device_loc']
    data = data.reset_index().set_index('id')
    # data_old = (pd.read_csv(os.path.join(data_dir, 'fiducial-desi4.dat'),
    #                         sep=r'\s+')
    #             .rename(columns={'#id': 'id'}).set_index('id'))

    data['x'] = [d['x'] for d in data[fid_centre_def]]
    data['y'] = [d['y'] for d in data[fid_centre_def]]
    data['z'] = [d['z'] for d in data[fid_centre_def]]
    # data['q_old'] = data_old['q']
    # data['s_old'] = data_old['s']
    data['q_kent'] = data_kent['Q']
    data['q_kent_wrapped'] = data['q_kent'] - (data['q_kent'] > 180) * 360
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
        data['q'], data['s'] = QS[0, :], QS[1, :]
        # data['dq_old_new'] = data['q_new'] - data['q_old']
        # data['ds_old_new'] = data['s_new'] - data['s_old']
        # data['dq_old'] = data['q_old'] - data['q_kent']
        # data['ds_old'] = data['s_old'] - data['s_kent']
        data['dq'] = data['q'] - data['q_kent_wrapped']
        data['ds'] = data['s'] - data['s_kent']
        # import pdb; pdb.set_trace()
        return np.sum(np.square(data[['dq', 'ds']].values))

    # minimise
    p0 = np.array([0, 0, 0, 0, 0, np.pi/5*(petal_loc-3)])
    # calculate_qst(p0)
    # data.to_csv('qst_test.csv')
    bounds = ((-10, 10), (-10, 10), (-10, 10),
              (-np.pi/2, np.pi//2), (-np.pi/2, np.pi/2), (-2*np.pi, 2*np.pi))
    solution = minimize(calculate_qst, p0, bounds=bounds, method='SLSQP',
                        options={'disp': True, 'maxiter': 10000})
    print(f'PTL{petal_id} Best-fit transformation params:\n{solution.x}')
    calculate_qst(solution.x)
    data.to_csv(f'qst_ptl{petal_id}_minimised.csv')
    return solution, data


if __name__ == '__main__':
    # petal_id = 11
    dfs = []
    rows = []
    for petal_id in range(3, 4):
        sol, data = qst_petal(petal_id)
        rows.append({'petal_id': petal_id,
                     'petal_loc': petal_locs[petal_id],
                     'offset_x': sol.x[0],
                     'offset_y': sol.x[1],
                     'offset_z': sol.x[2],
                     'rotation_x': sol.x[3],
                     'rotation_y': sol.x[4],
                     'rotation_z': sol.x[5],
                     'rotation_z_deg': np.degrees(sol.x[5]),
                     'least_square_residual': sol.fun})
        dfs.append(data)
    pd.DataFrame(rows).to_csv(os.path.join(
        data_dir, 'qst_focal_plate_alignment.csv'))
    pd.DataFrame(pd.concat(dfs, sort=False)).to_csv(os.path.join(
        data_dir, 'qst_fiducials_iter1.csv'))