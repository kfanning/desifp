# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 17:12:43 2020

@author: Duan Yutong (dyt@physics.bu.edu)
"""

import pandas as pd
from DOSlib.constants import ConstantsDB
from DOSlib.positioner_index import PositionerIndex
from posstate import PosState
from posmodel import PosModel

posid = 'M00001'
posintTP = (0, 135)
alt_calib = {'LENGTH_R1': 3.0,  # define calibration values here
             'LENGTH_R2': 3.0,
             'OFFSET_X': 0.0,
             'OFFSET_Y': 0.0,
             'OFFSET_T': 0.0,
             'OFFSET_P': 0.0}
pi = PositionerIndex()
# change tag to the time in history you want
# write a function to determine tag from DB given a python datetime object
group, tag, snapshot = 'focal_plane_metrology', get_tag(time), 'DESI'
df = pd.DataFrame.from_dict(ConstantsDB().get_constants(
    group=group, tag=tag, snapshot=snapshot)
    [group], orient='index')
alignments = {int(petal_loc): {'Tx': row['petal_offset_x'],
                               'Ty': row['petal_offset_y'],
                               'Tz': row['petal_offset_z'],
                               'alpha': row['petal_rot_1'],
                               'beta': row['petal_rot_2'],
                               'gamma': row['petal_rot_3']}
              for petal_loc, row in df.iterrows()}
state = PosState(unit_id=posid, device_type='pos')
petal_loc = pi.find_by_device_id(posid)['PETAL_LOC']
model = PosModel(state=state, petal_alignment=alignments[petal_loc])
trans = model.trans
trans.alt_override = True
trans.alt.update(alt_calib)
obsXY = trans.posintTP_to_obsXY(posintTP)
