# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 13:59:15 2019

@author: Duan Yutong (dyt@physics.bu.edu)
"""

import os
import numpy as np
import pandas as pd
from petaltransforms import PetalTransforms

fid_centre_def = 'center'
data_dir = r'/home/msdos/fp_temp_files/'
# data_dir = r'Downloads'
petal_locs = {2: 7, 3: 3, 4: 0, 5: 1, 6: 2, 7: 8, 8: 4, 9: 9, 10: 5, 11: 6}


def qst_petal(petal_id):
    data = pd.read_json(os.path.join(data_dir, f'petal{petal_id}.json'),
                        orient='index')
    petal_loc = petal_locs[petal_id]
    data['id'] = petal_loc*1000 + data['device_loc']
    data = data.reset_index().set_index('id')
    for coord in ['x', 'y', 'z']:
        data[coord] = [d[coord] for d in data[fid_centre_def]]
    # initialise petal transformation using nominal configuration, 36 deg * n
    trans = PetalTransforms(gamma=np.pi/5*(petal_loc-3))
    QST = trans.ptlXYZ_to_QST(data[['x', 'y', 'z']].values.T)
    data['q'], data['s'], data['t'] = QST[0, :], QST[1, :], QST[2, :]
    return data


if __name__ == '__main__':
    dfs = []
    for petal_id in petal_locs.keys():
        dfs.append(qst_petal(petal_id))
    pd.DataFrame(pd.concat(dfs, sort=False)).to_csv(os.path.join(
        data_dir, 'qst_fiducials-desi-dat.csv'))
