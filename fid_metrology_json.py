# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 13:21:35 2019

@author: Duan Yutong (dyt@physics.bu.edu)
"""
import os
import numpy as np
import pandas as pd
import simplejson as json

petal_id = 2
traveller_dir = r'Downloads'
traveller_fn = f'FPP Metrology Traveler - Petal{petal_id:02} - ZBF.xlsx'
path = os.path.join(traveller_dir, traveller_fn)
data = pd.read_excel(path, skiprows=17)
fifids = {11: 'P018', 75: 'P125', 150: 'P014', 239: 'P082', 321: 'P124',
          439: 'P104', 482: 'P010', 496: 'P066', 517: 'P023', 534: 'P096',
          541: 'P053', 542: 'P050'}

dump = {}
n_fifs = len(fifids)
for i in range(n_fifs):
    device_loc = data['Fiducial Locations'][1+i*5]
    fifid = fifids[device_loc]
    dump[fifid] = {'petal_id': int(petal_id),
                   'device_loc': device_loc}
    dump[fifid]['center'] = {
        'x': np.mean(data['As Built MeasuredZBF Locations'][2+i*5:6+i*5]),
        'y': np.mean(data['Unnamed: 7'][2+i*5:6+i*5]),
        'z': np.mean(data['Unnamed: 8'][2+i*5:6+i*5])}
with open(os.path.join(traveller_dir, f'petal{petal_id}.json'), 'w') as h:
    json.dump(dump, h, ensure_ascii=False, sort_keys=False, indent=4)
