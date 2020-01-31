# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 11:47:41 2020

@author: Duan Yutong (dyt@physics.bu.edu)

# set illumination level for each SpectroController
"""

import pandas as pd
from DOSlib.proxies import PetalMan, SimpleProxy

flashtime = 0.05  # second, 50 ms by default
data = {0: [0.77, 0.39, 1.00, 1.00],  # 0, 2, 3 are dim
        1: [1.00, 0.50, 0.25, 0.10],
        2: [0.50, 0.25, 0.30, 0.10],
        3: [0.65, 0.32, 0.50, 0.06],
        4: [0.80, 0.40, 0.30, 0.08],
        5: [1.00, 0.72, 0.34, 0.14],
        6: [1.00, 0.49, 0.40, 0.10],
        7: [1.00, 0.54, 1.00, 1.00],
        8: [0.49, 0.24, 0.53, 0.05],
        9: [1.00, 0.70, 0.43, 0.14]}
illum_cycles = pd.DataFrame.from_dict(data, orient='index',
                                      columns=[0.5, 1.0, 2.0, 5.0])
illum_cycles.index.name = 'PCID'
illum_cycles.columns.name = 'exptime/s'
fid_cycles = {0.5: 100, 1.0: 80, 2.0: 40, 5.0: 20}


def set_fiducials(enabled=True, exptime=2.0):
    ptlm = PetalMan()
    setting = fid_cycles[exptime] if enabled else 'off'
    retcode = ptlm.set_fiducials(setting=setting)
    print(f'PetalMan.set_fiducials {setting} returned: {retcode}')

def set_illuminators(enabled=True, exptime=2.0):
    if enabled:
        exptime = float(exptime)
        for pcid in illum_cycles.index:
            spectcon = SimpleProxy(f'SPECTCON{pcid}')
            spectcon._send_command('prepare_for_illumination')
            retcode = spectcon._send_command(
                'illuminator', action='on', flashtime=flashtime,
                duty_cycle=illum_cycles.loc[pcid, exptime])
            print(f'SPECTCON{pcid} on returned: {retcode}')
    else:
        specman = SimpleProxy('SPECMAN')
        retcode = specman._send_command('illuminate', action='off')
        print(f'SPECMAN off returned: {retcode}')

def rehome_adc():
    adc = SimpleProxy('ADC')
    retcode = adc._send_command('home', controllers=[1, 2])
    print(f'ADC.home returned: {retcode}')


if __name__ == '__main__':
    exptime = 2.0
    # exptime = float(input('FVC exposure time (s): '))
    action = input('Set fiducials (on/off, blank to skip): ')
    if action:
        assert action in {'on', 'off'}, 'Invalid input'
        enabled = True if action == 'on' else False
        set_fiducials(enabled=enabled, exptime=exptime)
    action = input('Set illuminators (on/off, blank to skip): ')
    if action:
        assert action in {'on', 'off'}, 'Invalid input'
        enabled = True if action == 'on' else False
        set_illuminators(enabled=enabled, exptime=exptime)
    action = input('Rehome ADC (y/n): ')
    if action == 'y':
        rehome_adc()
