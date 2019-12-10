# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 09:14:03 2019

@author: givoltage
"""
import pandas as pd
import numpy as np
path0 = r'Downloads\0.txt'
path1 = r'Downloads\1.txt'

data0 = pd.read_csv(path0, sep='\s+')
data1 = pd.read_csv(path1, sep='\s+')
data0['rot_deg'] = data0['ROT(DEG)'] + 36*(data0['PETAL_LOC']-3)
data0['rot_rad'] = np.radians(data0['rot_deg'])
data0['diff'] = data0['rot_rad'] - data1['petal_rot_3']
data0['diff mod 2pi'] = data0['diff'] % (2*np.pi)

data0['petal_rot_3'] = data1['petal_rot_3']
data0[['PETAL_LOC', 'ROT(DEG)', 'rot_rad', 'petal_rot_3', 'diff']].set_index('PETAL_LOC')
