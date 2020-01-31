# -*- coding: utf-8 -*-
"""
Created on Wed Jan  1 21:33:11 2020

@author: Duan Yutong

generator of json exposure script for FVC and Spotmatch stability tests
the resulting json file contains requests for all exptimes, N exposures each
run this json once for each telescope/dome configuration

cooling refers to:

    * Primary Mirror Thermal Control (PMTC):
        tl;dr:
            * only during the day with dome closed when telescope is not in use
            * mirror covers should be closed but could be open

        This draws dome air in through a filter (you can see the large hose on
        the SE side just under the railing on the C floor) into a heat
        exchanger located behind the south journal which valves building
        glycol to cool the air depending on the primary temperature setting
        and dew point.  The air is then blown up through the East elliptical
        tube into the space between the primary mirror and mirror covers.
        This is used only during the day when the telescope is not being used.
        One can hear the vibration from the large blower when it is on.

    * Primary Cell Vent Fans (PCVF):
        tl;dr:
            * night time dome closed only
            * mirror covers should be open

        These are the four squirrel cage fans on the center box section of the
        telescope at the NE, SE, SW, and NW locations.  These are controlled
        by the OA and are used to pull air from the space in front of the
        primary around the periphery of the primary mirror and out into the
        dome. The purpose of these is to avoid stagnant air cells which can
        form in the closed space within the center section in front of the
        primary and improve the image quality. As noted below, imaging
        experiments have not suggested any degradation of the image quality
        (and things should actually be better), but as far as I know, no high
        frequency measurements have ever been carried out. Whether these can
        induce vibration of the FVC is an open question. As noted, the OA can
        turn these on or off as desired.  They are used at night with the
        mirror covers open; the PMTC is always turned off when observing.

Other possible factors:
    * Dome Louvers:  While not properly fans, these are opened to flush air
    through the dome to minimize temperature differences between the telescope
    structure and the ambient air.  Depending on the wind speed and direction,
    it may be possible to set up resonances which could vibrate the telescope
    structure; this could be tested by opening/closing the louvers.

    * Dome Fans:  There are two large fans on the SSW side of the dome on the
    M floor level.  Prior to the installation of the louvers, this provided
    the only method of flushing air through the dome, and they were effectively
    useless in that respect (we carried out a test with neutral buoyancy
    balloons which demonstrated that turning these noisy fans on and off made
    no difference to the air flow in the dome).  I don't know if they are used
    anymore, but they do generate a bit of noise which could couple into
    telescope vibrations.

    * B29 Fan:  This may or may not be re-installed.  This is a huge fan which
    sits on the M floor to the East of the hatch and, as its name suggests, is
    extremely noisy.  It generates a vertical flow within the closed dome
    during the day to avoid thermal stratification of the air and is supposed
    to facilitate flushing of the dome air once the dome and louvers are opened
    for nighttime observing.  The noise makes any work within the dome during
    the day impossible, and this has typically been turned off when any daytime
    work is going on or the mirror covers are open (because of dust that could
    be stirred up).

"""

from itertools import product, chain
import simplejson as json

# exptimes = [0.5, 0.8, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]
exptimes = [0.5, 0.8, 1, 1.5, 2, 2.5, 3]  # keep under 3 to avoid saturation
conditions = {'dome': ['open', 'closed'],
              'angle': ['00', '25', '50'],
              'cooling': ['on', 'off']}
n_rep = 10  # repeat 10 times for each exptime
for cond in product(*conditions.values()):
    # condition is a tuple (dome, cooling, angle, ...)
    passthru = dict(zip(conditions.keys(), cond))
    fn = [f'{key} {val}' for key, val in zip(conditions.keys(), cond)]
    fn_ = [f'{key}_{val}' for key, val in zip(conditions.keys(), cond)]
    script = []
    for exptime in exptimes:
        script.append({'sequence': 'FVC', 'action': 'calibrate',
                       'exptime': exptime,
                       'program': f'fvc/sm stability: calibrate for '
                                  f'exptime {exptime}s, bias included',
                       'leave_fiducials': 'on', 'leave_illuminator': 'on',
                       'passthru': passthru})
        for i in range(n_rep):
            script.append(
                {'sequence': 'FVC', 'action': 'measure', 'exptime': exptime,
                 'program': f'fvc/sm stability: {", ".join(fn)} '
                            f'({i+1} of {n_rep})',
                 'fiducials': 'on', 'illuminator': 'on',
                 'leave_fiducials': 'on', 'leave_illuminator': 'on',
                 'passthru': passthru})
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
#             'fiducials': 'on', 'illuminator': 'on'}  # turn on for all exposures
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
