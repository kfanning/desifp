# -*- coding: utf-8 -*-
"""
Created on Tue Apr 9 19:05:20 2017

@author: Duan Yutong (dyt@physics.bu.edu)

Code on SVN:
    https://desi.lbl.gov/svn/code/focalplane/plate_layout/trunk/metrology_analysis/
Code on Github:
    https://github.com/givoltage/desifp
Results on google drive:
    https://drive.google.com/open?id=0B8mEgjxcgZeYVVFVZnYyM2Q3REE
Data on DocDB:
    DESI-3542: Focal plate survey and alignment data
    DESI-3543: Zeiss petal metrology data

"""

import os
from functools import reduce
# from multiprocessing import Pool
from concurrent.futures import ProcessPoolExecutor as Pool
import numpy as np
import pandas as pd
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages

metrology_dir = r'K:\Google Drive\DESI\model_drawings\DESI Focal Plate Assy and Integration\FP Structure\metrology\Zeiss Petal Metrology'
path_hole_table = r'K:\Google Drive\DESI\model_drawings\petal\DESI-0326-D_HoleTable.xlsx'
fig_save_dir = r'D:\plots'
petal_ids = ['01', '02', '04', '00', '03', '05', '06', '07', '08', '09',
             '10', '11']
paths = {'01': metrology_dir + r'\ptl01\DESI-0326-E_1_chr.txt',
         '02': metrology_dir + r'\ptl02\DESI-0326-E_2_chr.txt',
         '04': metrology_dir + r'\ptl04\DESI-0326-E_3_chr.txt',
         '00': metrology_dir + r'\ptl00\DESI-0326-E_5_chr.txt',
         '03': metrology_dir + r'\ptl03\DESI-0326-E_5_chr.txt',
         '05': metrology_dir + r'\ptl05\DESI-0326-E_6_chr.txt',
         '06': metrology_dir + r'\ptl06\DESI-0326-E_7_chr.txt',
         '07': metrology_dir + r'\ptl07\DESI-0326-E_8_chr.txt',
         '08': metrology_dir + r'\ptl08\DESI-0326-E_9_chr.txt',
         '09': metrology_dir + r'\ptl09\DESI-0326-E_10_chr.txt',
         '10': metrology_dir + r'\ptl10\DESI-0326-E_11_chr.txt',
         '11': metrology_dir + r'\ptl11\DESI-0326-E_12_chr.txt'}
make_plots = True
make_lookup_table = True

# %% function definitions


# rotation matrices, input is in radians
def Rx(angle):
    Rx = np.array([
                   [1.0,           0.0,            0.0],
                   [0.0,           np.cos(angle),  -np.sin(angle)],
                   [0.0,           np.sin(angle),  np.cos(angle)]
                  ])
    return Rx


def Ry(angle):
    Ry = np.array([
                   [np.cos(angle),  0.0,            np.sin(angle)],
                   [0.0,            1.0,            0.0],
                   [-np.sin(angle), 0.0,            np.cos(angle)]
                  ])
    return Ry


def Rz(angle):
    Rz = np.array([
                   [np.cos(angle), -np.sin(angle), 0.0],
                   [np.sin(angle), np.cos(angle),  0.0],
                   [0.0,           0.0,            1.0]
                  ])
    return Rz


def matmul(*args):
    return reduce(np.dot, [*args])


def total_residue(x, y):
    return np.sum(np.linalg.norm(x - y, axis=0))


def angles_to_unit_vector(theta_deg, phi_deg):
    # returns a unit vector from polar and azimuthal angles
    theta = np.radians(theta_deg)
    phi = np.radians(phi_deg)
    return np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])


def vector_to_angles(x, y, z):
    r = np.sqrt(np.square(x) + np.square(y) + np.square(z))
    theta = np.arccos(z/r) # radians
    phi = np.arctan(y/x)
    return [np.degrees(theta), np.degrees(phi)]


def petal_to_cs5(x, petal_location):
    
    '''
    given a petal location (designated by interger 0 to 9) and measurement 
    of a point in CS6, rotate the point by an integer number of 36 degrees
    to petal's local CS defined in the petal solid model.
    
    according to DESI-0742v3, the petal at location 3 shares the same 
    coordinate system with CS6 (X5 Y5 Z5)
    '''
    
    angle = 2*np.pi/10*(petal_location-3) # the point is rotated
    x_rot = matmul(Rz(angle), x)
    return x_rot


def throughput_tilt(tilt):
    # takes degree input
    throughput_tilt = -0.0133*np.power(tilt, 2) - 0.0175*tilt + 1.0 # takes deg
    return np.multiply(throughput_tilt, throughput_tilt>0)


def throughput_defocus(delta_f):
    delta_f_um = np.abs(delta_f) * 1000 # polyfit only defined for abs value
    # takes micron input
    throughput_defocus = (- 1.804e-14*np.power(delta_f_um, 5) # takes micron input
                          + 1.593e-11*np.power(delta_f_um, 4)
                          - 5.955e-10*np.power(delta_f_um, 3)
                          - 3.433e-6*np.power(delta_f_um, 2)
                          + 3.251e-7*delta_f_um
                          + 1.0)
    return np.multiply(throughput_defocus, throughput_defocus>0) # force positive valued


def par_list_gen(optm, half_ranges, counts):
    hr_fine = half_ranges[0]
    hr_medium = half_ranges[1]
    hr_coarse = half_ranges[2]
    array = np.concatenate((
                np.linspace(optm - hr_coarse, optm - hr_medium, counts[0]),
                np.linspace(optm - hr_medium, optm - hr_fine, counts[1]),
                np.linspace(optm - hr_fine, optm + hr_fine, counts[2]),
                np.linspace(optm + hr_fine, optm + hr_medium, counts[3]),
                np.linspace(optm + hr_medium, optm + hr_coarse, counts[4])))
    array.sort()
    return np.unique(array)

#%% analyse metrology data
# this function is going to be called by pool

def process_petal(petal_id):
    
    path = paths[petal_id]
    # petal_dir = os.path.dirname(path)
    df_txt = pd.read_table(path) # petal metrology data from txt
    df_hole_table = pd.read_excel(path_hole_table, skiprows=17)
    hole_ids = df_hole_table['hole_id']
    
    # transformation file
    df_transformations = pd.DataFrame(index = ['ABC', 'ZBF', 'SPT', 'TPT', '1DF', 'ACT'],
                                     columns = ['alpha', 'beta', 'gamma', 'Tx', 'Ty', 'Tz'])
    
    # initialise dataframe
    iterables = [['diameter', 'x', 'y', 'z', 'nutation', 'precession', 'r', 'tilt', 'defocus', 'throughput', 'dtb_x', 'dtb_y', 'dtb_z'],
                 ['ABC', 'ZBF', 'SPT', 'TPT', '1DF', 'ACT'],
                 hole_ids
                ]
    index = pd.MultiIndex.from_product(iterables, names=['hole feature', 'alignment', 'hole id'])
    df = pd.DataFrame(index=index, columns = df_txt.columns).sort_index()
    idx = pd.IndexSlice
    
    # ZBF alignment
    # fill in data from metrology report for 
    # ABC
    df.loc[idx['diameter', 'ABC']] = df_txt[df_txt['id'].str.contains('Sh5ZnC6 Diameter\(', na=False)].values
    df.loc[idx['x', 'ABC']] = df_txt[df_txt['id'].str.contains('Spotface Center X Position ABC Alignment\(', na=False)].values
    df.loc[idx['y', 'ABC']] = df_txt[df_txt['id'].str.contains('Spotface Center Y Position ABC Alignment\(', na=False)].values
    df.loc[idx['z', 'ABC']] = df_txt[df_txt['id'].str.contains('Spotface Center Z Position ABC Alignment\(', na=False)].values
    df.loc[idx['nutation', 'ABC']] = df_txt[df_txt['id'].str.contains('Nutation ABC Alignment\(', na=False)].values
    df.loc[idx['precession', 'ABC']] = df_txt[df_txt['id'].str.contains('Precession ABC Alignment\(', na=False)].values
    # ZBF
    df.loc[idx['diameter', 'ZBF']] = df_txt[df_txt['id'].str.contains('Sh5ZnC6 Diameter\(', na=False)].values
    df.loc[idx['x', 'ZBF']] = df_txt[df_txt['id'].str.contains('Spotface Center X Position Best Fit\(', na=False)].values
    df.loc[idx['y', 'ZBF']] = df_txt[df_txt['id'].str.contains('Spotface Center Y Position Best Fit\(', na=False)].values
    df.loc[idx['z', 'ZBF']] = df_txt[df_txt['id'].str.contains('Spotface Center Z Position Best Fit\(', na=False)].values
    df.loc[idx['nutation', 'ZBF']] = df_txt[df_txt['id'].str.contains('Nutation Best Fit\(', na=False)].values
    df.loc[idx['precession', 'ZBF']] = df_txt[df_txt['id'].str.contains('Precession Best Fit\(', na=False)].values
    # datum tooilng balls
    df.loc[('dtb_x', 'ZBF')][:3] = df_txt[df_txt['id'].str.contains('X Value_ToolingBall', na=False)].values
    df.loc[('dtb_y', 'ZBF')][:3] = df_txt[df_txt['id'].str.contains('Y Value_ToolingBall', na=False)].values
    df.loc[('dtb_z', 'ZBF')][:3] = df_txt[df_txt['id'].str.contains('Z Value_ToolingBall', na=False)].values

    # calculate r
    #ABC
    df.loc[idx['r', 'ABC'], 'nominal'] = np.sqrt(
                                np.square(df.loc[idx['x', 'ABC'], 'nominal'].values.astype(np.float64))
                                + np.square(df.loc[idx['y', 'ABC'], 'nominal'].values.astype(np.float64)))
    df.loc[idx['r', 'ABC'], 'actual'] = np.sqrt(
                                np.square(df.loc[idx['x', 'ABC'], 'actual'].values.astype(np.float64))
                                + np.square(df.loc[idx['y', 'ABC'], 'actual'].values.astype(np.float64)))
    #ZBF
    df.loc[idx['r', 'ZBF'], 'nominal'] = np.sqrt(
                                np.square(df.loc[idx['x', 'ZBF'], 'nominal'].values.astype(np.float64))
                                + np.square(df.loc[idx['y', 'ZBF'], 'nominal'].values.astype(np.float64)))
    df.loc[idx['r', 'ZBF'], 'actual'] = np.sqrt(
                                np.square(df.loc[idx['x', 'ZBF'], 'actual'].values.astype(np.float64))
                                + np.square(df.loc[idx['y', 'ZBF'], 'actual'].values.astype(np.float64)))
    
    # calculate tilt
    theta0 = df.loc['nutation', 'ZBF']['nominal'].values.astype(np.float64)
    phi0 = df.loc['precession', 'ZBF']['nominal'].values.astype(np.float64)
    theta = df.loc['nutation', 'ZBF']['actual'].values.astype(np.float64)
    phi = df.loc['precession', 'ZBF']['actual'].values.astype(np.float64)
    
    n0 = angles_to_unit_vector(theta0, phi0)
    n = angles_to_unit_vector(theta, phi)
    tilt = np.array([np.degrees(np.arccos(np.dot(n0[:, i], n[:, i]))) for i in range(514)])
    df.loc[idx['tilt', 'ZBF'], 'actual'] = tilt
    
    # calculate defocus
    delta_x = np.array([
                df.loc[idx['x', 'ZBF'], 'deviation'].values.astype(np.float64),
                df.loc[idx['y', 'ZBF'], 'deviation'].values.astype(np.float64),
                df.loc[idx['z', 'ZBF'], 'deviation'].values.astype(np.float64)
                ])
    delta_r = delta_x + 86.5 * (n-n0)
    delta_f = np.array([np.dot(delta_r[:, i], n0[:, i]) for i in range(514)])
    df.loc[idx['defocus', 'ZBF'], 'actual'] = delta_f
    df.loc['defocus', 'ABC']['lowertol'] = df.loc['z', 'ABC']['lowertol']
    df.loc['defocus', 'ABC']['uppertol'] = df.loc['z', 'ABC']['uppertol']
    df.loc[('tilt', 'ABC'), 'lowertol'] = 0
    df.loc[('tilt', 'ABC'), 'uppertol'] = 0.06
    df.loc[('throughput', 'ABC'), 'lowertol'] = 1-0.005 # 1 - throughput loss
    df.loc[('throughput', 'ABC'), 'uppertol'] = 1 # for throughput loss
    
    # calculate combined throughput    
    throughput = throughput_tilt(tilt) * throughput_defocus(delta_f)
    df.loc[idx['throughput', 'ZBF'], 'actual'] = throughput
    df.loc[idx['throughput', 'ZBF'], 'deviation'] = 1 - throughput

    #%% calculate transformation parameters of Zeiss' best-fit alignment

    x_abc = np.concatenate((df.loc[idx['x', 'ABC'],'actual'].values.reshape(-1, 1),
                            df.loc[idx['y', 'ABC'],'actual'].values.reshape(-1, 1),
                            df.loc[idx['z', 'ABC'],'actual'].values.reshape(-1, 1)),
                            axis=1).T.astype(np.float64)
    x_zbf = np.concatenate((df.loc[idx['x', 'ZBF'],'actual'].values.reshape(-1, 1),
                            df.loc[idx['y', 'ZBF'],'actual'].values.reshape(-1, 1),
                            df.loc[idx['z', 'ZBF'],'actual'].values.reshape(-1, 1)),
                            axis=1).T.astype(np.float64)
    x_nom = np.concatenate((df.loc[idx['x', 'ZBF'],'nominal'].values.reshape(-1, 1),
                            df.loc[idx['y', 'ZBF'],'nominal'].values.reshape(-1, 1),
                            df.loc[idx['z', 'ZBF'],'nominal'].values.reshape(-1, 1)),
                            axis=1).T.astype(np.float64)
    theta_abc = np.radians(df.loc[idx['nutation', 'ABC'],'actual'].values.astype(np.float64))
    phi_abc = np.radians(df.loc[idx['precession', 'ABC'],'actual'].values.astype(np.float64))
    x0_abc = x_abc - np.array([np.sin(theta_abc)*np.cos(phi_abc),
                               np.sin(theta_abc)*np.sin(phi_abc),
                               np.cos(theta_abc)])  # another point along axis
    
    def total_residue_1(parameters):
        # this is the function to be minimised
        alpha = parameters[0]
        beta = parameters[1]
        gamma = parameters[2]
        T = parameters[3:]
        R = matmul(Rz(gamma), Ry(beta), Rx(alpha)) # yaw-pitch-roll system
        x_abc_rot = np.empty(x_abc.shape) # rotated matrix from ABC
        # rotate each column of x_abc and fill x_abc_rot
        for j in range(x_abc.shape[1]):
            x_abc_rot[:, j] = matmul(R, x_abc[:, j]) + T
        
        return total_residue(x_abc_rot, x_zbf)
    
    # minimisation routine
    p0 = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    # alpha rotation about z-axis is limited by the 0.6 mm gap between petals
    bounds = ((-0.0026, 0.0026), (-np.pi/2, np.pi/2), (-np.pi/2, np.pi/2), (-10, 10), (-10, 10), (-10, 10))
    solution_zbf = minimize(total_residue_1, p0, 
                          bounds = bounds, 
                          method='SLSQP', 
                          options={'disp': True})
#    # mod out 2pi
#    for k in range(3):
#        parameters[k] = parameters[k] % (2*np.pi)
    print('PTL'+petal_id+' ZBF Transformation recovered as \n'
          +'X Rotation (Roll) : {}° \n'.format(np.degrees(solution_zbf.x[0]))
          +'Y Rotation (Pitch): {}° \n'.format(np.degrees(solution_zbf.x[1]))
          +'Z Rotation (Yaw)  : {}° \n'.format(np.degrees(solution_zbf.x[2]))
          +'Translation: {} \n'.format(solution_zbf.x[3:])
          +'With least square: {} \n'.format(solution_zbf.fun))
    
    # save transformation parameters
    df_transformations.loc['ZBF'] = solution_zbf.x # in radians
    
    # Calculate R and T
    parameters = solution_zbf.x
    alpha = parameters[0]
    beta = parameters[1]
    gamma = parameters[2]
    T = parameters[3:]
    R = matmul(Rz(gamma), Ry(beta), Rx(alpha)) # yaw-pitch-roll system
    
    # datum tooilng balls
    for j in range(3):
        
        dtb_pos_zbf = np.array([df.loc[('dtb_x', 'ZBF'), 'actual'][j],
                                df.loc[('dtb_y', 'ZBF'), 'actual'][j],
                                df.loc[('dtb_z', 'ZBF'), 'actual'][j]])
        dtb_pos_abc = matmul(np.linalg.inv(R), dtb_pos_zbf - T)
        df.loc[('dtb_x', 'ABC'), 'actual'][j] = dtb_pos_abc[0]
        df.loc[('dtb_y', 'ABC'), 'actual'][j] = dtb_pos_abc[1]
        df.loc[('dtb_z', 'ABC'), 'actual'][j] = dtb_pos_abc[2]
    
    #%% 'SPT' alignment
    # geometric fit, using spotface centres only
    
    df.loc[('diameter', 'SPT')] = df.loc[('diameter', 'ABC')].values
    
    def total_residue_2(parameters):
        # this is the function to be minimised
        alpha = parameters[0] # radians
        beta = parameters[1]
        gamma = parameters[2]
        T = parameters[3:]
        R = matmul(Rz(gamma), Ry(beta), Rx(alpha)) # yaw-pitch-roll system
        x_abc_rot = np.empty(x_abc.shape) # rotated matrix from ABC
        # rotate each column of x_abc and fill x_abc_rot
        for j in range(x_abc.shape[1]):
            x_abc_rot[:, j] = matmul(R, x_abc[:, j]) + T
        return total_residue(x_abc_rot, x_nom)

    solution_spt = minimize(total_residue_2, p0, 
                            bounds = bounds, 
                            method='SLSQP', 
                            options={'disp': True})
    
    print('PTL'+petal_id+' SPT transformation found as \n'
          +'X Rotation (Roll) : {}° \n'.format(np.degrees(solution_spt.x[0]))
          +'Y Rotation (Pitch): {}° \n'.format(np.degrees(solution_spt.x[1]))
          +'Z Rotation (Yaw)  : {}° \n'.format(np.degrees(solution_spt.x[2]))
          +'Translation: {} \n'.format(solution_spt.x[3:])
          +'With least square: {} \n'.format(solution_spt.fun))
    
    # save transformation parameters
    df_transformations.loc['SPT'] = solution_spt.x
    
    # write results with these parameters to dataframe
    parameters = solution_spt.x
    alpha = parameters[0] # radians
    beta = parameters[1]
    gamma = parameters[2]
    T = parameters[3:]
    R = matmul(Rz(gamma), Ry(beta), Rx(alpha)) # yaw-pitch-roll system
    x_abc_rot = np.empty(x_abc.shape) # rotated matrix from ABC
    x0_abc_rot = np.empty(x_abc.shape)

    # rotate each column of x_abc and fill x_abc_rot
    for j in range(x_abc.shape[1]):
        x_abc_rot[:, j] = matmul(R, x_abc[:, j]) + T
        x0_abc_rot[:, j] = matmul(R, x0_abc[:, j]) + T
    
    n_rot = x_abc_rot - x0_abc_rot
    [theta_rot, phi_rot] = vector_to_angles(n_rot[0, :], n_rot[1, :], n_rot[2, :])
                
    df.loc[idx['x', 'SPT'], 'actual'] = x_abc_rot[0, :]
    df.loc[idx['y', 'SPT'], 'actual'] = x_abc_rot[1, :]
    df.loc[idx['z', 'SPT'], 'actual'] = x_abc_rot[2, :]
    df.loc[idx['nutation', 'SPT'], 'actual'] = theta_rot
    df.loc[idx['precession', 'SPT'], 'actual'] = phi_rot
    df.loc[idx['x', 'SPT'], 'nominal'] = df.loc[idx['x', 'ABC'], 'nominal'].values
    df.loc[idx['y', 'SPT'], 'nominal'] = df.loc[idx['y', 'ABC'], 'nominal'].values
    df.loc[idx['z', 'SPT'], 'nominal'] = df.loc[idx['z', 'ABC'], 'nominal'].values
    df.loc[idx['nutation', 'SPT'], 'nominal'] = df.loc[idx['nutation', 'ABC'], 'nominal'].values
    df.loc[idx['precession', 'SPT'], 'nominal'] = df.loc[idx['precession', 'ABC'], 'nominal'].values
    df.loc[idx['x', 'SPT'], 'deviation'] = df.loc[idx['x', 'SPT'], 'actual'].values - df.loc[idx['x', 'ABC'], 'nominal'].values
    df.loc[idx['y', 'SPT'], 'deviation'] = df.loc[idx['y', 'SPT'], 'actual'].values - df.loc[idx['y', 'ABC'], 'nominal'].values
    df.loc[idx['z', 'SPT'], 'deviation'] = df.loc[idx['z', 'SPT'], 'actual'].values - df.loc[idx['z', 'ABC'], 'nominal'].values
    df.loc[idx['nutation', 'SPT'], 'deviation'] = df.loc[idx['nutation', 'SPT'], 'actual'].values - df.loc[idx['nutation', 'ABC'], 'nominal'].values
    df.loc[idx['precession', 'SPT'], 'deviation'] = df.loc[idx['precession', 'SPT'], 'actual'].values - df.loc[idx['precession', 'ABC'], 'nominal'].values
    
    # calculate r 
    df.loc[idx['r', 'SPT'], 'actual'] = np.sqrt(
                            np.square(df.loc[idx['x', 'SPT'], 'actual'].values.astype(np.float64))
                            + np.square(df.loc[idx['y', 'SPT'], 'actual'].values.astype(np.float64)))
    
    # datum tooilng balls
    for j in range(3):
        x = np.array([df.loc[('dtb_x', 'ABC'), 'actual'][j],
                      df.loc[('dtb_y', 'ABC'), 'actual'][j],
                      df.loc[('dtb_z', 'ABC'), 'actual'][j]])
        xp = matmul(R, x) + T
        df.loc[('dtb_x', 'SPT'), 'actual'][j] = xp[0]
        df.loc[('dtb_y', 'SPT'), 'actual'][j] = xp[1]
        df.loc[('dtb_z', 'SPT'), 'actual'][j] = xp[2]
    
    # calculate tilt
    theta = df.loc['nutation', 'SPT']['actual'].values.astype(np.float64)
    phi = df.loc['precession', 'SPT']['actual'].values.astype(np.float64)
    
    n = angles_to_unit_vector(theta, phi)
    tilt = np.array([np.degrees(np.arccos(np.dot(n0[:, i], n[:, i]))) for i in range(514)])
    df.loc[idx['tilt', 'SPT'], 'actual'] = tilt
    
    # calculate defocus
    delta_x = np.array([
                df.loc[idx['x', 'SPT'], 'deviation'].values.astype(np.float64),
                df.loc[idx['y', 'SPT'], 'deviation'].values.astype(np.float64),
                df.loc[idx['z', 'SPT'], 'deviation'].values.astype(np.float64)
                ])
    delta_r = delta_x + 86.5 * (n-n0)
    delta_f = np.array([np.dot(delta_r[:, i], n0[:, i]) for i in range(514)])
    df.loc[idx['defocus', 'SPT'], 'actual'] = delta_f
    
    # calculate combined throughput
    throughput = throughput_tilt(tilt) * throughput_defocus(delta_f)
    df.loc[idx['throughput', 'SPT'], 'actual'] = throughput
    df.loc[idx['throughput', 'SPT'], 'deviation'] = 1 - throughput
    
    #%% 'TPT' alignment, throughput-based fitting
    
    df.loc[('diameter', 'TPT')] = df.loc[('diameter', 'ABC')].values
    
    def throughput_loss_min(parameters):
        # this is the function to be minimised
        alpha = parameters[0] # radians
        beta = parameters[1]
        gamma = parameters[2]
        T = parameters[3:]
        R = matmul(Rz(gamma), Ry(beta), Rx(alpha)) # yaw-pitch-roll system
        x_abc_rot = np.empty(x_abc.shape) # rotated matrix from ABC
        x0_abc_rot = np.empty(x_abc.shape)
        # rotate each column of x_abc and fill x_abc_rot
        for k in range(x_abc.shape[1]):
            x_abc_rot[:, k] = matmul(R, x_abc[:, k]) + T
            x0_abc_rot[:, k] = matmul(R, x0_abc[:, k]) + T
        n_rot = x_abc_rot - x0_abc_rot # not necessarily of unit length
        [theta_rot, phi_rot] = vector_to_angles(n_rot[0, :], n_rot[1, :], n_rot[2, :])
        n_rot = angles_to_unit_vector(theta_rot, phi_rot) # renormalise to unit length
        # calculate defocus
        delta_x = x_abc_rot - x_nom
        delta_r = delta_x + 86.5 * (n_rot-n0)
        # calculate tilt
        tilt = np.array([np.degrees(np.arccos(np.dot(n0[:, i], n_rot[:, i]))) for i in range(514)])
        delta_f = np.array([np.dot(delta_r[:, i], n0[:, i]) for i in range(514)])
        throughput = throughput_tilt(tilt) * throughput_defocus(delta_f)
        
        return np.mean(1-throughput)

    solution_tpt = minimize(throughput_loss_min, p0,
                          bounds = bounds, 
                          method = 'SLSQP', 
                          options={'disp': True,
                                   'maxiter': 10000})
    
    print('PTL'+petal_id+' TPT transformation found as \n'
          +'X Rotation (Roll) : {}° \n'.format(np.degrees(solution_tpt.x[0]))
          +'Y Rotation (Pitch): {}° \n'.format(np.degrees(solution_tpt.x[1]))
          +'Z Rotation (Yaw)  : {}° \n'.format(np.degrees(solution_tpt.x[2]))
          +'Translation: {} \n'.format(solution_tpt.x[3:])
          +'With least mean throughput loss: {} \n'.format(solution_tpt.fun))
    
    # save transformation parameters
    df_transformations.loc['TPT'] = solution_tpt.x
    
    # write results with these parameters to dataframe
    parameters = solution_tpt.x
    alpha = parameters[0] # radians
    beta = parameters[1]
    gamma = parameters[2]
    T = parameters[3:]
    R = matmul(Rz(gamma), Ry(beta), Rx(alpha)) # yaw-pitch-roll system
    x_abc_rot = np.empty(x_abc.shape) # rotated matrix from ABC
    x0_abc_rot = np.empty(x_abc.shape)
    
    # rotate each column of x_abc and fill x_abc_rot
    for j in range(x_abc.shape[1]):
        x_abc_rot[:, j] = matmul(R, x_abc[:, j]) + T
        x0_abc_rot[:, j] = matmul(R, x0_abc[:, j]) + T
        
    n_rot = x_abc_rot - x0_abc_rot
    [theta_rot, phi_rot] = vector_to_angles(n_rot[0, :], n_rot[1, :], n_rot[2, :])
    
    df.loc[('x', 'TPT'), 'actual'] = x_abc_rot[0, :]
    df.loc[('y', 'TPT'), 'actual'] = x_abc_rot[1, :]
    df.loc[('z', 'TPT'), 'actual'] = x_abc_rot[2, :]
    df.loc[('nutation', 'TPT'), 'actual'] = theta_rot
    df.loc[('precession', 'TPT'), 'actual'] = phi_rot
    df.loc[('x', 'TPT'), 'nominal'] = df.loc[('x', 'ABC'), 'nominal'].values
    df.loc[('y', 'TPT'), 'nominal'] = df.loc[('y', 'ABC'), 'nominal'].values
    df.loc[('z', 'TPT'), 'nominal'] = df.loc[('z', 'ABC'), 'nominal'].values
    df.loc[('nutation', 'TPT'), 'nominal'] = df.loc[('nutation', 'ABC'), 'nominal'].values
    df.loc[('precession', 'TPT'), 'nominal'] = df.loc[('precession', 'ABC'), 'nominal'].values
    df.loc[('x', 'TPT'), 'deviation'] = df.loc[('x', 'TPT'), 'actual'].values - df.loc[('x', 'ABC'), 'nominal'].values
    df.loc[('y', 'TPT'), 'deviation'] = df.loc[('y', 'TPT'), 'actual'].values - df.loc[('y', 'ABC'), 'nominal'].values
    df.loc[('z', 'TPT'), 'deviation'] = df.loc[('z', 'TPT'), 'actual'].values - df.loc[('z', 'ABC'), 'nominal'].values
    df.loc[('nutation', 'TPT'), 'deviation'] = df.loc[('nutation', 'TPT'), 'actual'].values - df.loc[('nutation', 'ABC'), 'nominal'].values
    df.loc[('precession', 'TPT'), 'deviation'] = df.loc[('precession', 'TPT'), 'actual'].values - df.loc[('precession', 'ABC'), 'nominal'].values
    
    # datum tooilng balls
    for j in range(3):
        x = np.array([df.loc[('dtb_x', 'ABC'), 'actual'][j],
                      df.loc[('dtb_y', 'ABC'), 'actual'][j],
                      df.loc[('dtb_z', 'ABC'), 'actual'][j]])
        xp = matmul(R, x) + T
        df.loc[('dtb_x', 'TPT'), 'actual'][j] = xp[0]
        df.loc[('dtb_y', 'TPT'), 'actual'][j] = xp[1]
        df.loc[('dtb_z', 'TPT'), 'actual'][j] = xp[2]
    
    # calculate r 
    df.loc[idx['r', 'TPT'], 'actual'] = np.sqrt(
                            np.square(df.loc[idx['x', 'TPT'], 'actual'].values.astype(np.float64))
                            + np.square(df.loc[idx['y', 'TPT'], 'actual'].values.astype(np.float64)))
    
    # calculate tilt
    theta = df.loc['nutation', 'TPT']['actual'].values.astype(np.float64)
    phi = df.loc['precession', 'TPT']['actual'].values.astype(np.float64)
    
    n = angles_to_unit_vector(theta, phi)
    tilt = np.array([np.degrees(np.arccos(np.dot(n0[:, i], n[:, i]))) for i in range(514)])
    df.loc[idx['tilt', 'TPT'], 'actual'] = tilt
    
    # calculate defocus
    delta_x = np.array([
                df.loc[idx['x', 'TPT'], 'deviation'].values.astype(np.float64),
                df.loc[idx['y', 'TPT'], 'deviation'].values.astype(np.float64),
                df.loc[idx['z', 'TPT'], 'deviation'].values.astype(np.float64)
                ])
    delta_r = delta_x + 86.5 * (n-n0)
    delta_f = np.array([np.dot(delta_r[:, i], n0[:, i]) for i in range(514)])
    df.loc[idx['defocus', 'TPT'], 'actual'] = delta_f
    
    # calculate combined throughput
    throughput = throughput_tilt(tilt) * throughput_defocus(delta_f)
    df.loc[idx['throughput', 'TPT'], 'actual'] = throughput
    df.loc[idx['throughput', 'TPT'], 'deviation'] = 1 - throughput
    
    #%% ACT alignment, actual tooling ball positions measured on CMM
    
    df.loc[('x', 'ACT'), 'nominal'] = df.loc[('x', 'ABC'), 'nominal'].values
    df.loc[('y', 'ACT'), 'nominal'] = df.loc[('y', 'ABC'), 'nominal'].values
    df.loc[('z', 'ACT'), 'nominal'] = df.loc[('z', 'ABC'), 'nominal'].values
    df.loc[('nutation', 'ACT'), 'nominal'] = df.loc[('nutation', 'ABC'), 'nominal'].values
    df.loc[('precession', 'ACT'), 'nominal'] = df.loc[('precession', 'ABC'), 'nominal'].values
    
    #%% all plots
    if make_plots:
        
        figtitle_prefix = 'Petal ' + str(petal_id)
        figtitles = {'diameter': 'Cylinder Diameter Deviation',
                     'x': 'Spotface Centre X Deviation',
                     'y': 'Spotface Centre Y Deviation',
                     'z': 'Spotface Centre Z Deviation',
                     'nutation': 'Nutation Deviation',
                     'precession': 'Precession Deviation',
                     'tilt': 'Cylinder Axial Tilt',
                     'defocus': 'Spotface Centre Defocus',
                     'throughput': 'Throughput Loss'}
        axtitles = {'diameter': r'$\delta D/\mathrm{mm}$',
                    'x': r'$\delta x/\mathrm{mm}$',
                    'y': r'$\delta y/\mathrm{mm}$',
                    'z': r'$\delta z/\mathrm{mm}$',
                    'nutation': r'$\delta \theta/\degree$',
                    'precession': r'$\delta \varphi/\degree$',
                    'tilt': r'$\delta/\degree$',
                    'defocus': r'$\delta f/\mathrm{mm}$',
                    'throughput': r'$1-\eta$'}
        colour_range = {'diameter': [0.008, 0.018],
                        'x': [-0.03, 0.03],
                        'y': [-0.03, 0.03],
                        'z': [-0.03, 0.03],
                        'nutation': [-0.028, 0.028],
                        'precession': [-0.028, 0.028],
                        'tilt': [0, 0.05],
                        'defocus': [-0.03, 0.03],
                        'throughput': [0, 0.005]}
        tol_lower = {'diameter': df.loc['diameter', 'ABC']['lowertol'],
                     'x': df.loc['x', 'ABC']['lowertol'],
                     'y': df.loc['y', 'ABC']['lowertol'],
                     'z': df.loc['z', 'ABC']['lowertol'],
                     'nutation': df.loc['nutation', 'ABC']['lowertol'],
                     'precession': df.loc['precession', 'ABC']['lowertol'],
                     'tilt': df.loc['tilt', 'ABC']['lowertol'],
                     'defocus': df.loc['defocus', 'ABC']['lowertol'],
                     'throughput': 1-df.loc['throughput', 'ABC']['lowertol']}
        tol_upper = {'diameter': df.loc['diameter', 'ABC']['uppertol'],
                     'x': df.loc['x', 'ABC']['uppertol'],
                     'y': df.loc['y', 'ABC']['uppertol'],
                     'z': df.loc['z', 'ABC']['uppertol'],
                     'nutation': df.loc['nutation', 'ABC']['uppertol'],
                     'precession': df.loc['precession', 'ABC']['uppertol'],
                     'tilt': df.loc['tilt', 'ABC']['uppertol'],
                     'defocus': df.loc['defocus', 'ABC']['uppertol'],
                     'throughput': 1-df.loc['throughput', 'ABC']['uppertol']}
        units = {'diameter': ' mm',
                 'x': ' mm',
                 'y': ' mm',
                 'z': ' mm',
                 'nutation': r'$\degree$',
                 'precession': r'$\degree$',
                 'tilt': r'$\degree$',
                 'defocus': ' mm',
                 'throughput': ''}
    
        for alignment in ['ZBF', 'SPT', 'TPT']:
            
            x = df.loc['x', alignment]['nominal']
            y = df.loc['y', alignment]['nominal']
            r = df.loc['r', alignment]['actual']
            colours = {'diameter': df.loc['diameter', alignment]['deviation'],
                       'x': df.loc['x', alignment]['deviation'],
                       'y': df.loc['y', alignment]['deviation'],
                       'z': df.loc['z', alignment]['deviation'],
                       'nutation': df.loc['nutation', alignment]['deviation'],
                       'precession': df.loc['precession', alignment]['deviation'],
                       'tilt': df.loc['tilt', alignment]['actual'],
                       'defocus': df.loc['defocus', alignment]['actual'],
                       'throughput': df.loc['throughput', alignment]['deviation']}
            
            for feature in ['diameter', 'x', 'y', 'z', 'nutation', 'precession', 'tilt', 'defocus', 'throughput']:
                
                fig = plt.figure(figsize = (18, 6))
                gs = gridspec.GridSpec(1, 2, width_ratios=[3, 2])
                fig.suptitle(figtitle_prefix + ' ' + figtitles[feature] + ' (' + alignment + ' Alignment)')
                rms = np.sqrt(np.mean(np.square(colours[feature])))
                mean = np.mean(colours[feature])
                sd = np.std(colours[feature])
                textstr = (r'RMS={0:.7f}' + units[feature] + '\n'
                           + r'$\mu={1:.7f}$' + units[feature] + '\n'
                           + r'$\sigma={2:.7f}$'
                           + units[feature]).format(rms, mean, sd)
                textbbox = {'boxstyle': 'square',
                            'facecolor': 'white',
                            'alpha': 0.5}
                ax0 = plt.subplot(gs[0])
                plot0 = ax0.scatter(x, y, marker='o', lw = 4, 
                                    c=colours[feature], cmap='plasma', 
                                    vmin = colour_range[feature][0], vmax = colour_range[feature][1])
                fig.colorbar(plot0, ax=ax0)
                ax0.axis('equal')
                ax0.set_title(axtitles[feature])
                ax0.set_xlabel('x/mm')
                ax0.set_ylabel('y/mm')
                ax0.text(50, 150, textstr, fontsize=12, linespacing=1.5, bbox=textbbox, ha='left', va='bottom')
                ax1 = plt.subplot(gs[1])
                ax1.plot(r, colours[feature], 'b.', r, tol_upper[feature], 'r--', r, tol_lower[feature], 'r--', )
                ax1.set_xlabel('r/mm')
                ax1.set_ylabel(axtitles[feature])
                fig.tight_layout()
                # save figure
                figname = str(petal_ids.index(petal_id)+1).zfill(2) \
                          + '-ptl_' + petal_id + '-' + alignment.lower() \
                          + '-' + feature + '.pdf'
                # fig.savefig(os.path.join(fig_save_dir, figname), dpi=600)
                pp = PdfPages(os.path.join(fig_save_dir, figname))
                pp.savefig(fig, dpi=1000)
                pp.close()
                plt.close('all')
    else:
        pass

    # export dataframes
    df.to_csv(os.path.join(fig_save_dir, 
                           str(petal_ids.index(petal_id)+1).zfill(2) \
                           + '-ptl_' + petal_id + '-df_data.csv'))
    df.to_pickle(os.path.join(fig_save_dir, 
                           str(petal_ids.index(petal_id)+1).zfill(2) \
                           + '-ptl_' + petal_id + '-df_data.pickle'),
                 compression='gzip')
    df_transformations.to_csv(os.path.join(fig_save_dir, 
                                          str(petal_ids.index(petal_id)+1).zfill(2) \
                                          + '-ptl_' + petal_id + '-df_transformations.csv'))
    df_transformations.to_pickle(os.path.join(fig_save_dir, 
                                          str(petal_ids.index(petal_id)+1).zfill(2) \
                                          + '-ptl_' + petal_id + '-df_transformations.pickle'),
                                compression='gzip')
    
    #%% ICS conformity per DESI-2850
    # everything in petal local CS
    
    device_loc = np.append(hole_ids.values, [543, 544, 545])
    x_meas = np.append(df.loc['x','ZBF']['actual'].values,
                       df.loc['dtb_x','ZBF']['actual'][0:3].values)
    y_meas = np.append(df.loc['y','ZBF']['actual'].values,
                       df.loc['dtb_y','ZBF']['actual'][0:3].values)
    z_meas = np.append(df.loc['z','ZBF']['actual'].values,
                       df.loc['dtb_z','ZBF']['actual'][0:3].values)
    x_meas_proj = np.append(df.loc['x','ZBF']['actual'].values + 86.5*n[0,:], 
                            [np.nan]*3)
    y_meas_proj = np.append(df.loc['y','ZBF']['actual'].values + 86.5*n[1,:],
                            [np.nan]*3)
    z_meas_proj = np.append(df.loc['z','ZBF']['actual'].values + 86.5*n[2,:],
                            [np.nan]*3)
    
    df_ics = pd.DataFrame(
                {'petal_id': petal_id,
                 'device_loc': device_loc,
                 'x_meas': x_meas,
                 'y_meas': y_meas,
                 'z_meas': z_meas,
                 'x_meas_proj': x_meas_proj,
                 'y_meas_proj': y_meas_proj,
                 'z_meas_proj': z_meas_proj
                 })
    df_ics.to_csv(os.path.join(fig_save_dir, 
                           str(petal_ids.index(petal_id)+1).zfill(2) \
                           + '-ptl_' + petal_id + '-petal_metrology_ics.csv'),
                  index = False, index_label = False, 
                  columns = ['petal_id', 'device_loc', 
                             'x_meas', 'y_meas', 'z_meas',
                             'x_meas_proj', 'y_meas_proj', 'z_meas_proj'])
    df_ics.to_pickle(os.path.join(fig_save_dir, 
                           str(petal_ids.index(petal_id)+1).zfill(2) \
                           + '-ptl_' + petal_id + '-petal_metrology_ics.pickle'),
                     compression = 'gzip')
            
    #%% DTB coordinates for all possible mounting positions
    iterables = [np.arange(10), ['dtb_0', 'dtb_1', 'dtb_2']]
    index = pd.MultiIndex.from_product(iterables, names = ['ring mounting location', 'dtb id'])
    df_dtb = pd.DataFrame(index = index, columns = ['x', 'y', 'z']).sort_index()
    for i in np.arange(10):
        for j, dtb_id in enumerate(['dtb_0', 'dtb_1', 'dtb_2']):
            x = np.array([df.loc['dtb_x','ZBF']['actual'][j],
                          df.loc['dtb_y','ZBF']['actual'][j],
                          df.loc['dtb_z','ZBF']['actual'][j]])
            df_dtb.loc[i, dtb_id] = petal_to_cs5(x, i)
        # iterate over 10 possible mounting locations on the ring
    df_dtb.to_csv(os.path.join(fig_save_dir, 
                           str(petal_ids.index(petal_id)+1).zfill(2) \
                           + '-ptl_' + petal_id + '-df_dtb_optimal.csv'))
    df_dtb.to_pickle(os.path.join(fig_save_dir, 
                           str(petal_ids.index(petal_id)+1).zfill(2) \
                           + '-ptl_' + petal_id + '-df_dtb_optimal.pickle'),
                 compression='gzip')
    
    #%% pre-calculated throughput table in petal local ABC coordinate system
    # while varying around ZBF parameters as initial values
    
    def tpt_eval(alpha, beta, gamma, Tx, Ty, Tz):
        # given transforamtion parameters w.r.t. ABC, return throughput
        T = np.array([Tx, Ty, Tz])
        R = matmul(Rz(gamma), Ry(beta), Rx(alpha)) # yaw-pitch-roll system
        
        # calculate tooling ball positions
        dtb_pos = np.empty((3,3))
        for i in range(3):
            x = np.array([df.loc[('dtb_x', 'ABC'), 'actual'][i],
                          df.loc[('dtb_y', 'ABC'), 'actual'][i],
                          df.loc[('dtb_z', 'ABC'), 'actual'][i]])
            xp = matmul(R, x) + T
            dtb_pos[i, :] = xp
        # calculate throughput
        x_abc_rot = np.empty(x_abc.shape) # rotated matrix from ABC
        x0_abc_rot = np.empty(x_abc.shape)
        for j in range(x_abc.shape[1]):
            x_abc_rot[:, j] = matmul(R, x_abc[:, j]) + T
            x0_abc_rot[:, j] = matmul(R, x0_abc[:, j]) + T
        n_rot = x_abc_rot - x0_abc_rot
        [theta_rot, phi_rot] = vector_to_angles(n_rot[0, :], n_rot[1, :], n_rot[2, :])
        n_rot = angles_to_unit_vector(theta_rot, phi_rot)
        tilt = np.array([np.degrees(np.arccos(np.dot(n0[:, i], n_rot[:, i]))) for i in range(514)])
        delta_x = x_abc_rot - np.array(
                [df.loc[('x', 'ABC'), 'nominal'].values.astype(np.float64),
                 df.loc[('y', 'ABC'), 'nominal'].values.astype(np.float64),
                 df.loc[('z', 'ABC'), 'nominal'].values.astype(np.float64)])
        delta_r = delta_x + 86.5 * (n-n0)
        delta_f = np.array([np.dot(delta_r[:, i], n0[:, i]) for i in range(514)])
        throughput = throughput_tilt(tilt) * throughput_defocus(delta_f)
    
        return (dtb_pos[0, 0], dtb_pos[0, 1], dtb_pos[0, 2],
                dtb_pos[1, 0], dtb_pos[1, 1], dtb_pos[1, 2],
                dtb_pos[2, 0], dtb_pos[2, 1], dtb_pos[2, 2],
                np.mean(throughput), np.std(throughput), 
                np.sqrt(np.mean(np.square(throughput))),
                np.mean(1-throughput), np.std(1-throughput),
                np.sqrt(np.mean(np.square(1-throughput))))
    
    if make_lookup_table:
        
        tpt_eval_v = np.vectorize(tpt_eval)
        
    #    parameters_list = zip(np.arange(alpha0 - 1e-4, alpha0 + 1e-4, 5e-5),
    #            np.arange(beta0 - 5e-4, beta0 + 5e-4, 2e-4),
    #            np.arange(gamma0 - 1e-3, gamma0 + 1e-3, 5e-4),
    #            np.arange(Tx0 - 0.1, Tx0 + 0.1, 0.05),
    #            np.arange(Ty0 - 0.1, Ty0 + 0.1, 0.05),
    #            np.arange(Tz0 - 0.1, Tz0 + 0.1, 0.05))
    #    pool_tpt_lookup = Pool()
    #    results = pool_tpt_lookup.starmap(tpt_eval, parameters_list)
    #    rows = [result[0] for result in results]
    #    df_tpt_lookup = pd.concat(rows)
        
        alpha0 = df_transformations.loc['ZBF']['alpha']
        beta0 = df_transformations.loc['ZBF']['beta']
        gamma0 = df_transformations.loc['ZBF']['gamma']
        Tx0 = df_transformations.loc['ZBF']['Tx']
        Ty0 = df_transformations.loc['ZBF']['Ty']
        Tz0 = df_transformations.loc['ZBF']['Tz']
    
        alpha = par_list_gen(alpha0, [1e-4, 4e-4, 8e-4], [1, 1, 1, 1, 1])
        beta = par_list_gen(beta0, [2e-4, 6e-4, 12e-4], [1, 4, 10, 4, 1]) #10
        gamma = par_list_gen(gamma0, [3e-4, 8e-4, 8e-4], [0, 1, 1, 1, 0])
        Tx = par_list_gen(Tx0, [0.1, 0.3, 0.3], [0, 1, 1, 1, 0])
        Ty = par_list_gen(Ty0, [0.1, 0.3, 0.3], [0, 1, 1, 1, 0])
        Tz = par_list_gen(Tz0, [0.05, 0.15, 0.3], [1, 4, 12, 4, 1]) #12
        # note that meshgrid axis=0 correspond to 2nd argument, beta, 
        # axis=1 corresonds to 1st argument, alpha. weird bug
        alpha_mg, beta_mg, gamma_mg, Tx_mg, Ty_mg, Tz_mg = np.meshgrid(
                alpha, beta, gamma, Tx, Ty, Tz)
        print('PTL'+petal_id+' Total number of entries in parameter space for' \
              + ' PTL{} throughput evaluation: {} \n'.format(petal_id, alpha_mg.size))
        results = tpt_eval_v(alpha_mg, beta_mg, gamma_mg, Tx_mg, Ty_mg, Tz_mg)
        # switch beta and alpha when taking gradient to match dimensions
        gradient = np.gradient(results[9], beta, alpha, gamma, Tx, Ty, Tz)
        df_tpt_lookup = pd.DataFrame({'alpha': alpha_mg.flatten(),
                                      'beta':  beta_mg.flatten(),
                                      'gamma': gamma_mg.flatten(),
                                      'Tx': Tx_mg.flatten(),
                                      'Ty': Ty_mg.flatten(),
                                      'Tz': Tz_mg.flatten(),
                                      'dtb_0_x': results[0].flatten(),
                                      'dtb_0_y': results[1].flatten(),
                                      'dtb_0_z': results[2].flatten(),
                                      'dtb_1_x': results[3].flatten(),
                                      'dtb_1_y': results[4].flatten(),
                                      'dtb_1_z': results[5].flatten(),
                                      'dtb_2_x': results[6].flatten(),
                                      'dtb_2_y': results[7].flatten(),
                                      'dtb_2_z': results[8].flatten(),
                                      'throughput_mean': results[9].flatten(),
                                      'throughput_std': results[10].flatten(),
                                      'throughput_rms': results[11].flatten(),
                                      'throughput_loss_mean': results[12].flatten(),
                                      'throughput_loss_std': results[13].flatten(),
                                      'throughput_loss_rms': results[14].flatten(),
                                      'grad_alpha': gradient[1].flatten(),
                                      'grad_beta': gradient[0].flatten(),
                                      'grad_gamma': gradient[2].flatten(),
                                      'grad_x': gradient[3].flatten(),
                                      'grad_y': gradient[4].flatten(),
                                      'grad_z': gradient[5].flatten()})
        
        df_tpt_lookup.to_csv(os.path.join(fig_save_dir, 
                               str(petal_ids.index(petal_id)+1).zfill(2) \
                               + '-ptl_' + petal_id + '-throughput_lookup.csv'))
        df_tpt_lookup.to_pickle(os.path.join(fig_save_dir, 
                               str(petal_ids.index(petal_id)+1).zfill(2) \
                               + '-ptl_' + petal_id + '-throughput_lookup.pickle'),
                     compression='gzip')
        print('PTL'+petal_id+' TPT lookup table saved')
    else:
        pass

if __name__ == '__main__':
    
    # determine quality of production petal by analysing Zeiss metrology data
    pool = Pool() # create a multiprocessing Pool
    pool.map(process_petal, petal_ids) # process Zeiss metrology data
    