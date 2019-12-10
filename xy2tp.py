# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 19:14:46 2019

@author: givoltage
"""

import sys
import math
import posconstants as pc


r1 = 3.02794845598864
r2 = 3.04335035678846
r = [r1, r2]
targetable_range_T = [-193.277142730862, 190.277142730862]
targetable_range_P = [2.2696165167000117, 182.0]
ranges = [targetable_range_T, targetable_range_P]


def tp2xy(tp, r):
    """Converts TP angles into XY cartesian coordinates, where arm lengths
    associated with angles theta and phi are respectively r[1] and r[2].
    INPUTS:  tp ... [theta,phi], unit degrees
              r ... [central arm length, eccentric arm length]
    OUTPUT:  xy ... [x,y]
    """
    t = math.radians(tp[0])
    t_plus_p = t + math.radians(tp[1])
    x = r[0] * math.cos(t) + r[1] * math.cos(t_plus_p)
    y = r[0] * math.sin(t) + r[1] * math.sin(t_plus_p)
    return x, y


def xy2tp(xy, r, ranges):
    """Converts XY cartesian coordinates into TP angles, where arm lengths
     associated with angles theta and phi are respectively r[1] and r[2].

    INPUTS:   xy ... [x,y]
               r ... [central arm length, eccentric arm length]
          ranges ... [[min(theta), max(theta)], [min(phi), max(phi)]]

    OUTPUTS:  tp ... [theta,phi], unit degrees
     unreachable ... boolean, True if the requested xy cannot be reached
                     by any tp

    In cases where unreachable == True, the returned tp value will be a
    closest possible approach to the unreachable point requested at xy.
    """
    # within this much xy error allowance, adjust theta toward center
    # of its range
    theta_centralizing_err_tol = 1e-4
    # number of points to try when attempting to centralize theta
    n_theta_centralizing_iters = 3
    # slight contraction to avoid numeric divide-by-zero type of errors
    numeric_contraction = sys.float_info.epsilon*10
    x, y, r1, r2 = xy[0], xy[1], r[0], r[1]
    unreachable = False
    # adjust targets within reachable annulus
    hypot = (x**2.0 + y**2.0)**0.5
    angle = math.atan2(y, x)
    outer = r[0] + r[1]
    inner = abs(r[0] - r[1])
    if hypot > outer or hypot < inner:
        unreachable = True
    inner += numeric_contraction
    outer -= numeric_contraction
    HYPOT = hypot
    if hypot >= outer:
        HYPOT = outer
    elif hypot <= inner:
        HYPOT = inner
    X = HYPOT*math.cos(angle)
    Y = HYPOT*math.sin(angle)
    # transfrom from cartesian XY to angles TP
    arccos_arg = (X**2.0 + Y**2.0 - (r1**2.0 + r2**2.0)) / (2.0 * r1 * r2)
    # deal with slight numeric errors where arccos_arg comes back
    # like -1.0000000000000002
    arccos_arg = max(arccos_arg, -1.0)
    # deal with slight numeric errors where arccos_arg comes back
    # like +1.0000000000000002
    arccos_arg = min(arccos_arg, +1.0)
    P = math.acos(arccos_arg)
    T = angle - math.atan2(r2*math.sin(P), r1 + r2*math.cos(P))
    TP = [math.degrees(T), math.degrees(P)]
    # wrap angles into travel ranges
    for i in [0, 1]:
        range_min, range_max = min(ranges[i]), max(ranges[i])
        if TP[i] < range_min:
            # try +360 phase wrap
            TP[i] += math.floor((range_max - TP[i])/360.0)*360.0
            if TP[i] < range_min:
                # print(f'TP {i}, {TP[i]} < range_min: {range_min}'
                #       f'thus unreachable\n{ranges[i]}')
                TP[i] = range_min
                unreachable = True
        elif TP[i] > range_max:
            # try -360 phase wrap
            TP[i] -= math.floor((TP[i] - range_min)/360.0)*360.0
            if TP[i] > range_max:
                # print(f'TP {i}, {TP[i]} > range_max: {range_max}'
                #       f'thus unreachable\n{ranges[i]}')
                TP[i] = range_max
                unreachable = True
    # centralize theta
    T_ctr = (ranges[0][0] + ranges[0][1])/2.0
    T_options = pc.linspace(TP[0], T_ctr, n_theta_centralizing_iters)
    for T_try in T_options:
        xy_try = tp2xy([T_try, TP[1]], r)  # change to ref class staticmethod
        x_err = xy_try[0] - X
        y_err = xy_try[1] - Y
        vector_err = (x_err**2.0 + y_err**2.0)**0.5
        if vector_err <= theta_centralizing_err_tol:
            TP[0] = T_try
            break
    return tuple(TP), unreachable


xy2tp((-0.015401900799819845, 3.7270292731249426e-16), r, ranges)
