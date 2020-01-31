# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 17:41:59 2020

@author: Duan Yutong (dyt@physics.bu.edu)

  - Telescope at zenith
  - Nighttime, dome open
  - 104 fiducials only, with manual duty cycle control
  - Exposures: 10x 0.5, 1, 2, 5s exposures
  - Take one set with "Primary Cell Vent Fans" on, another set with them off.
  - These should be taken in succession to minimize the impact of other
    variables.

may have to use a dummy PCID in pecs_local config if empty list doesn't work

"""
import pandas as pd
from configobj import ConfigObj
import posconstants as pc
from pecs import PECS
from fptestdata import FPTestData

exptimes = [0.5, 0.8, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]
duties = {0.5: 100, 0.8: 100, 1: 100, 1.5: 80, 2: 60, 2.5: 50,
          3: 40, 3.5: 35, 4: 30, 4.5: 25, 5: 20}
n_fid = 104
test = PECS(interactive=True)
cfg = ConfigObj()
cfg['pcids'] = test.pcids
test.data = FPTestData('FVC stability', test_cfg=cfg)
test.logger = test.data.logger  # broadcast to all petals
test.data.t_i = pc.now()
test.exp_setup()
exp_exptimes, mea_exptimes = [], []  # one giant dataframe for each exptime
test.fvc.make_targets(num_spots=4*n_fid)  # make FVC targets using # pinholes
for exptime in exptimes:
    test.print(f'Setting exptime to {exptime} s')
    duty = duties[exptime]
    test.ptlm.set_fiducials(setting=duty)  # set fiducial duty cycle
    test.fvc.calibrate_bias(dark_flag=1)  # take FVC dark bias for new exptime
    test.fvc.calibrate_image()  # might retry here depending on what's returned
    n_rep = 100 if exptime == 0.5 else 10
    exp_frames, mea_frames, = [], []  # list of dataframes to be concatenated
    for i in range(n_rep):
        test.print(f'Taking exposure {i+1} of {n_rep}')
        exppos, meapos, _, _ = test.fvc_measure()
        exp_frames.append(exppos)
        mea_frames.append(meapos)
    exp_exptimes.append(pd.concat(exp_frames, axis=1, keys=range(n_rep)))
    mea_exptimes.append(pd.concat(exp_frames, axis=1, keys=range(n_rep)))
test.data.t_f = pc.now()
test.exppos = pd.concat(exp_exptimes, axis=1, keys=exptimes,
                        names=['exptime', 'expnum', 'field'])
test.meapos = pd.concat(mea_exptimes, axis=1, keys=exptimes,
                        names=['exptime', 'expnum', 'field'])
test.fvc_collect(destination=test.data.dir)
test.data.export_data_logs()
test.data.dump_as_one_pickle()
test.data.make_archive()
