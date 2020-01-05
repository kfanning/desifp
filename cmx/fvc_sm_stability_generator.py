# -*- coding: utf-8 -*-
"""
Created on Wed Jan  1 21:33:11 2020

@author: Duan Yutong

generator of json exposure script for FVC and Spotmatch stability tests
the resulting json file contains requests for all exptimes, N exposures each
run this json once for each telescope/dome configuration
"""
from itertools import product, chain
import simplejson as json

exptimes = [0.5, 0.8, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]
conditions = {'dome': ['open', 'closed'],
              'angle': ['00', '25', '50'],
              'cooling': ['on', 'off']}
n_rep = 10  # repeat 10 times for each exptime
for cond in product(*conditions.values()):
    # condition is a tuple (dome, cooling, angle, ...)
    fn = [f'{key} {val}' for key, val in zip(conditions.keys(), cond)]
    fn_ = [f'{key}_{val}' for key, val in zip(conditions.keys(), cond)]
    script = []
    for exptime in exptimes:
        script.append({'sequence': 'FVC', 'action': 'calibrate',
                       'exptime': exptime,
                       'program': f'fvc/sm stability: calibrate for '
                                  f'exptime {exptime}s, bias included',
                       'leave_fiducials': 'on', 'leave_illuminator': 'on'})
        for i in range(n_rep):
            script.append(
                {'sequence': 'FVC', 'action': 'measure', 'exptime': exptime,
                 'program': f'fvc/sm stability: {", ".join(fn)} '
                            f'({i+1} of {n_rep})',
                 'fiducial': 'on', 'illuminator': 'on',
                 'leave_fiducials': 'on', 'leave_illuminator': 'on'})
    with open(f'{"-".join(fn_)}-{n_rep}x.json', 'w') as h:
        json.dump(script, h, ensure_ascii=False, indent=4)


# conditions = {'dome': ['open', 'closed'],
#               'angle': ['00', '25', '50'],
#               'cooling': ['on', 'off'],
#               'exptime': [0.5, 0.8, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]}
# n_rep = 10  # repeat 10 times for each exptime
# for cond in product(*conditions.values()):
#     script = []
#     exptime = cond[3]
#     # condition is a tuple (dome, cooling, angle, ...)
#     fn = [f'{key} {val}' for key, val in zip(conditions.keys(), cond)]
#     fn_ = [f'{key}_{val}' for key, val in zip(conditions.keys(), cond)]
#     req = {'sequence': 'FVC', 'action': 'measure', 'exptime': exptime,
#             'program': f'fvc/sm stability: {", ".join(fn)}',
#             'fiducial': 'on', 'illuminator': 'on'}  # turn on for all exposures
#     for i in range(n_rep):
#         if i == 0:  # first image for new exposure time, calibrate first
#             script.append({'sequence': 'FVC', 'action': 'calibrate',
#                             'exptime': exptime,
#                             'program': f'fvc/sm stability: calibrate for '
#                                       f'{exptime}s after manual bias/dark',
#                             'leave_fiducials': 'on', 'leave_illuminator': 'on'})
#             req_i = copy(req)
#             req_i.update({'leave_fiducials': 'on', 'leave_illuminator': 'on'})
#             req_i['program'] = f'{req_i["program"]}, expnum {i+1} of {n_rep}'
#             script.append(req_i)
#         elif i == n_rep - 1:  # last img, turn off illuminator
#             req_f = copy(req)
#             req_f.update({'leave_fiducials': 'off', 'leave_illuminator': 'off'})
#             req_f['program'] = f'{req_f["program"]}, expnum {i+1} of {n_rep}'
#             script.append(req_f)
#         else:
#             req_m = copy(req)
#             req_m.update({'leave_fiducials': 'on', 'leave_illuminator': 'on'})
#             req_m['program'] = f'{req_m["program"]}, expnum {i+1} of {n_rep}'
#             script.append(req_m)
#     with open(f'{"-".join(fn_)}-{n_rep}x.json', 'w') as h:
#         json.dump(script, h, ensure_ascii=False, indent=4)
