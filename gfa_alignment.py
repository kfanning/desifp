# -*- coding: utf-8 -*-
"""
Created on Mon Apr  2 23:59:13 2018

@author: givoltage
"""

import numpy as np
from functools import reduce

petal_location_lookup = {'00': 2,  # map between petal location and petal ID
                         '01': 7,
                         '02': 7,
                         '03': 8,
                         '04': 0,
                         '05': 1,
                         '06': 2,
                         '07': 3,
                         '08': 4,
                         '09': 9,
                         '10': 5,
                         }
gfa_datums = {'00': np.array([[401.163, -35.0627, -98.9602],
                              [367.319, -30.0371, -95.7755],
                              [367.3465, -7.0439, -95.6652]]),
              '01': np.array([[-401.1384, 35.0217, -98.9542],
                              [-367.2983, 29.9563, -95.7617],
                              [-367.3518, 6.9623, -95.6577]]),
              '02': np.array([[-401.0465, 34.9637, -98.9393],
                              [-367.1945, 29.9708, -95.7387],
                              [-367.2043, 6.9751, -95.6338]]),
              '03': np.array([[-345.0852, -207.5725, -98.9417],
                              [-314.7758, -191.6971, -95.7616],
                              [-301.2542, -210.2968, -95.6453]]),
              '04': np.array([[90.5371, -392.2961, -98.9412],
                              [84.8286, -358.555, -95.7562],
                              [106.7119, -351.4668, -95.6377]]),
              '05': np.array([[303.9361, -264.1682, -98.9488],
                              [279.5309, -240.1882, -95.7563],
                              [293.0876, -221.6166, -95.6456]]),
              '06': np.array([[401.1716, -35.0441, -98.9527],
                              [367.3277, -29.9975, -95.7706],
                              [367.37, -7.0021, -95.6439]]),
              '07': np.array([[345.1574, 207.3516, -98.9326],
                              [314.7978, 191.5788, -95.7478],
                              [301.3364, 210.2213, -95.6454]]),
              '08': np.array([[157.3129, 370.6949, -98.9565],
                              [142.0578, 340.0719, -95.76],
                              [120.2016, 347.2139, -95.655]]),
              '09': np.array([[-157.2776, -370.7094, -98.9563],
                              [-142.0624, -340.0582, -95.7675],
                              [-120.1944, -347.1709, -95.6584]]),
              '10': np.array([[-90.5619, 392.3265, -98.9571],
                              [-84.93, 358.5781, -95.7761],
                              [-106.8184, 351.5298, -95.6584]])
              }


# %%
def matmul(*args):
    return reduce(np.dot, [*args])


def Rx(angle):
    Rx = np.array([
                   [1.0,           0.0,            0.0],
                   [0.0,           np.cos(angle),  np.sin(angle)],
                   [0.0,           -np.sin(angle), np.cos(angle)]
                  ])
    return Rx


def Ry(angle):
    Ry = np.array([
                   [np.cos(angle), 0.0,           -np.sin(angle)],
                   [0.0,           1.0,           0.0],
                   [np.sin(angle), 0.0,           np.cos(angle)]
                  ])
    return Ry


def Rz(angle):
    Rz = np.array([
                   [np.cos(angle),  np.sin(angle), 0.0],
                   [-np.sin(angle), np.cos(angle), 0.0],
                   [0.0,            0.0,           1.0]
                  ])
    return Rz


def Rxyz(alpha, beta, gamma):

    Rxyz = matmul(Rz(gamma), Ry(beta), Rx(alpha))  # yaw-pitch-roll system
    return Rxyz


def cs5_to_petal(x, petal_location):

    '''
    given a petal location (designated by interger 0 to 9) and measurement
    of a point in CS6, rotate the point by an integer number of 36 degrees
    to petal's local CS defined in the petal solid model.

    according to DESI-0742v3, the petal at location 3 shares the same
    coordinate system with CS6 (X5 Y5 Z5)
    '''

    theta = 2*np.pi/10*(petal_location-3)
    x_rot = matmul(Rz(theta), x)
    return x_rot


# %%
# GFA nominal datums positions at location 3
x_0 = np.array([[345.109, 207.422, -98.940],
                [314.792, 191.557, -95.743],
                [301.273, 210.164, -95.631]])
n = np.cross(x_0[0, :] - x_0[1, :], x_0[0, :] - x_0[2, :])
n = n/np.linalg.norm(n)
# cols: deltan1, deltan2, deltan3, rms delta n, residual sum of squares
delta_n = np.zeros((12, 5))
for j, petal_id in enumerate(petal_location_lookup.keys()):
    petal_location = petal_location_lookup[petal_id]
    x_cs5 = gfa_datums[petal_id]  # measured GFA datums positions in CS5
    x = np.zeros((3, 3))
    delta_r = []
    for i in range(3):
        x[i, :] = cs5_to_petal(x_cs5[i, :], petal_location)  # in local CS
        delta_x = x[i, :] - x_0[i, :]
        delta_n[int(petal_id), i] = np.dot(delta_x, n)
        delta_n[int(petal_id), 4] += np.square(np.linalg.norm(delta_x))
        delta_r.append(np.linalg.norm(delta_x[:2]))
    print('PTL {}, R err: {}, \n N err: {}'
          .format(petal_id, delta_r, delta_n[int(petal_id), :]))

delta_n[:, 3] = np.sqrt(np.mean(np.square(delta_n[:, :3]), axis=1))
delta_n[:, 4] = np.sqrt(delta_n[:, 4]/3)
print(delta_n)

