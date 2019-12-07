# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 16:15:02 2019

@author: Duan Yutong (dyt@physics.bu.edu)

Take FIF metrology data, FIF install info, and an FVC image, generate
FIF metrology data from the optical measurement.
Find magnification, translation, rotation angle which takes petal CS to
FVC pixel space
    x_fvc = M * (R(gamma) @ x_petal + T)
as well as the rotation for each FIF, assuming the measured cap position
lies exactly at the positioner hole centre at the focal surface

Assumptions:
    FIF is perfectly straight, all 4 pinholes lie in the same plane
    perpendicular to the chief-ray direction as defined in DESI-0530.
    The central axis of FIF completely overlaps with this chief-ray direction,
    and FIF can rotate by an arbitrary amount about this axis. This gives
    13 parameters from 24 points.

    Then holding these 13 parameters fixed, allow in-plane shifts of pinholes
    possibly due to FIF being not straight or bending. Now fit another 12
    parameters. Two iterations yield a total of 25 free parameters.
"""

import os
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from astropy.io import fits
import simplejson as json
from petal_metrology import Ry, Rz, Rxyz, total_residue
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family': 'serif',
                     'mathtext.fontset': 'cm'})


# fifid, device_loc, approx fvcX, fvcY
fif_info = {'P078': [439,  800, 2235],
            'P059': [482,  860, 2800],
            'P067': [517,  980, 3480],
            'P056': [321, 2120, 3140],
            'P002': [496, 1760, 4380],
            'P073': [534, 1300, 4270]}
data_dir = r'C:\Users\givoltage\Google Drive\DESI\model_drawings'
path_fif_met = os.path.join(
    data_dir, 'Production_DESI_Fiducial_Measurements2.xlsx')
path_petal_met = os.path.join(data_dir, 'Petal_Metrology.csv')
path_fvc_img = os.path.join(data_dir, 'fvc.20190717150817.fits')
path_fif_fvc = os.path.join(data_dir, 'fvc.20190717150817.pos')
fif_met = pd.read_excel(path_fif_met, skiprows=4)
petal_met = pd.read_csv(path_petal_met)
fif_fvc = np.loadtxt(path_fif_fvc)[:, :2]
img = fits.open(path_fvc_img)[0]
ydim, xdim = img.shape


def fvc_flip(xdata, fvc_yflip=-1):
    '''  # -1 because FVC does flip image about y-axis'''
    return fvc_yflip * xdata + (fvc_yflip < 0) * xdim


def check_colinear(positions):
    '''input 4 x 2 array, return point that's not colinear and others'''
    for i, pos in enumerate(positions):
        otherpos = np.vstack([positions[:i], positions[i+1:]])  # three others
        matrix = np.vstack([[1, 1, 1], otherpos.T])
        area = np.abs(np.linalg.det(matrix))  # area of triangle
        if 0 < area < 0.2:
            # area of 3 linearly independent points is about 1 or 2
            # if 3 colienar, area close to 0
            return positions[i], otherpos


def check_farthest(reference, positions):
    '''return the point in positions that's farthest to reference'''
    i = np.argmax(np.linalg.norm(positions - reference, axis=1))
    return positions[i], np.vstack([positions[:i], positions[i+1:]])


def identify_pinholes(positions):
    '''input 4 x 2 array, return ordered array according to pinhole numbers'''
    ordered = np.zeros(positions.shape)
    pt1, otherpos = check_colinear(positions)
    ordered[0] = pt1
    pt2, otherpos = check_farthest(pt1, otherpos)
    ordered[1] = pt2
    ordered[2:4, :] = np.vstack(check_farthest(pt2, otherpos))
    return ordered


def parse_fif_measurements(fif_fvc):
    '''parse xy positions (N x 2) from source extraction and associate 4
    points with each FIF using approximate position with 50 px radius'''
    fif_data_dict = {}
    for fifid, info in fif_info.items():
        pos_approx = np.array(info[1:])
        # look for points close to centre
        mask = np.linalg.norm(fif_fvc - pos_approx, axis=1) < 50
        fif_data_dict[fifid] = identify_pinholes(fif_fvc[mask])
    return fif_data_dict


def read_fif_metrology(fifid):
    '''read FIF local metrology data (2 x 5) in shifted local centre frame'''
    fifid_short = f'P{int(fifid[1:])}'
    idx = fif_met[fif_met['Fiducial Number'] == fifid_short].index
    x_fid = np.zeros((2, 5))  # in fid CS
    for i in range(5):
        x_fid[:, i] = fif_met.loc[idx+i, ['X (mm)', 'Y (mm)']].values
    # move origin to as-measured cap centre (point 5)
    x_fid = x_fid - x_fid[:, 4].reshape(2, 1)  # all 5 pinholes in 5 columns
    return x_fid  # in fiducial local CS shifted to cap centre


def Rzyz(device_loc):
    '''rotation matrix for transformation from fid/pos local CS to petal CS
    using Z(spin)Y(nutation)Z(precession) convention
    From DESI-0387 sheet 6:
        Step 1: rotate around Z by precession (spherical polar phi)
        Step 2: rotate around Y by nutation (spherical polar theta)
        Step 3: unrotate around z by spin
    '''
    theta = np.radians(petal_met.loc[device_loc, 'nutat'])
    phi = np.radians(petal_met.loc[device_loc, 'preces'])
    spin = np.radians(petal_met.loc[device_loc, 'spin'])
    return Rz(-spin) @ Ry(theta) @ Rz(phi)


def transform_model(params, inplane_shifts=None):
    '''inplane_shifts, if supplied, should be 2 x N_fif'''
    if inplane_shifts is None:
        inplane_shifts = np.zeros((3, len(fif_info)))
    else:  # add z=0 to 3rd dim, shape is now 3 x N_fif
        inplane_shifts = np.vstack([inplane_shifts.reshape(2, -1),
                                    np.zeros(inplane_shifts.shape[1])])
    alpha, beta, gamma = np.radians(params[:3])  # input angles are in degrees
    Tx, Ty, Tz, M = params[3:7]  # float PTL and FVC params
    angles = np.radians(params[7:])  # float FIF rotation params in radians
    # create model for FIF pinhole positions, with float params above
    fif_petal = {}  # in petal CS
    for i, fifid in enumerate(fif_info.keys()):
        device_loc = fif_info[fifid][0]  # look up nominal centre in petal CS
        # as-measured pinhole positions in local CS, z is 0 for all 5 positions
        x_fid = np.vstack([read_fif_metrology(fifid), np.zeros(5)])
        # rotate by arbitrary angle about local z axis
        x_fid = Rz(angles[i]) @ x_fid + inplane_shifts[:, i].reshape(3, 1)
        T = (petal_met.loc[device_loc, ['X', 'Y', 'Z']]  # cap centre offset
             .values.astype(np.float64).reshape(3, 1))
        # take pinhole positions from local CS to petal CS, still in 3D
        fif_petal[fifid] = (Rzyz(device_loc) @ x_fid + T).astype(np.float64)
    model = np.hstack([fif_petal[fifid][:, :4] for fifid in fif_info.keys()])
    # transform model to match observed data in FVC pixel space
    T = np.array([Tx, Ty, Tz]).reshape(3, 1)
    modelprime = M * (Rxyz(alpha, beta, gamma) @ model + T)
    return modelprime, fif_petal  # modelprime is FVC CS, model is petal CS


def reduced_chisq(deviations, covariance):
    '''use chisq if we have multiple exposures to establish covariance'''
    return deviations @ np.linalg.inv(covariance) @ deviations.T


def residue_ptl_fif(params, data):  # input angles are in degrees
    '''13 DoF OLS estimator accouting for petal rotatation, translation,
    FIF rotation, and FVC magnification (somewhat degenerate with Tz)'''
    assert len(params) == 7 + len(fif_info.keys())
    modelprime, _ = transform_model(params)
    # treat as two independent var, delta x (1st row) and delta y (2nd row)
    rss = total_residue(data, modelprime[:2, :])  # only 1 exposure, use OLS
    print(f'Residue func eval:  = {rss}')
    return rss


def residue_inplane_shifts(shifts, params, data):
    '''12 DoF OLS estimator accounting for local in-plane shifts of pinholes'''
    assert shifts.size == 2 * len(fif_info.keys())  # 2 x N describing 2D
    modelprime, _ = transform_model(params,
                                    inplane_shifts=np.reshape(shifts, (2, -1)))
    # treat as two independent var, delta x (1st row) and delta y (2nd row)
    rss = total_residue(data, modelprime[:2, :])  # only 1 exposure, use OLS
    print(f'Residue func eval: {rss}')
    return rss


def plot_solutions(solution1, solution2):
    '''plot 13 DoF and 12 DoF fits to show goodness and in-plane corrections'''
    text = (
        f'Best-fit transformation parameters ({len(solution1.x)} DoF):'
        f'\n3-component petal rotations: ' +
        (r'$(\alpha, \beta, \gamma) = '
         r'({: 8.3f}\degree, {: 8.3f}\degree, {: 8.3f}\degree)$').format(
            solution1.x[0], solution1.x[1], solution1.x[2]) +
        f'\n3-component petal translation: ' +
        (r'$\vec{{T}} = ({: 8.3f}\mathrm{{\,px}}, '
         r'{: 8.3f}\mathrm{{\,px}}, {: 8.3f}\mathrm{{\,px}})$').format(
         solution1.x[3], solution1.x[4], solution1.x[5]) +
        f'\n1 magnification: 'r'$M = {: 8.3f} \mathrm{{\,px/mm}}$'
        .format(solution1.x[6]) +
        f'\n{len(fif_info)} FIF additional local spins:')
    for i, fifid in enumerate(fif_info.keys()):
        text += '\n    'r'$\varphi(\mathrm{{{0}}}) = {1: 8.3f}\degree$'.format(
            fifid, solution1.x[7+i])
    text += ('\nLeast square sum of residues: '
             r'${: 8.3f} \mathrm{{\,px}}^2$'.format(solution1.fun) +
             f'\nWith additional {len(solution2.x)} DoF '
             'for in-plane shifts of pinholes: '
             r'${: 8.3f} \mathrm{{\,px}}^2$'.format(solution2.fun))
    modelprime1, fif_ptl1 = transform_model(solution1.x)
    modelprime2, fif_ptl2 = transform_model(
        solution1.x, inplane_shifts=np.reshape(solution2.x, (2, -1)))
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.plot(modelprime1[0], modelprime1[1], 'o', fillstyle='none',
            ms=3, mew=0.2, color='C3', label='Transformed nominal model')
    ax.plot(modelprime2[0], modelprime2[1], 'D', fillstyle='none',
            ms=2, mew=0.2, color='b', label='Model with in-plane shifts')
    ax.plot(fif_data[0], fif_data[1], 'k+', ms=2.5, mew=0.2,
            label='FVC measurements')
    for fifid, info in fif_info.items():
        ax.text(fvc_flip(info[1], fvc_yflip=-1)-66, info[2]-70,
                f'{fifid} ({info[0]})', fontsize=9)
    ax.text(0.18, 0.05, text, transform=ax.transAxes,  # family='monospace',
            fontsize=10, linespacing=2,
            horizontalalignment='left', verticalalignment='bottom',
            bbox={'boxstyle': 'round', 'alpha': 0.8,
                  'facecolor': 'white', 'edgecolor': 'lightgrey'})
    ax.set_xlabel(r'$x_\mathrm{FVC}/\mathrm{px}$')
    ax.set_ylabel(r'$y_\mathrm{FVC}/\mathrm{px}$')
    ax.legend(loc='upper right', fontsize=10, markerscale=3)
    fig.savefig('fif_ptl_fvc_metrology.pdf')
    return modelprime2, fif_ptl2


if __name__ == '__main__':
    fif_data_dict = parse_fif_measurements(fif_fvc)  # fvc pinhole positions
    fif_data = np.hstack([fif_data_dict[fifid].T for fifid in fif_info.keys()])
    # unflip FVC image
    fif_data[0] = fvc_flip(fif_data[0], fvc_yflip=-1)
    iv1 = np.array([-2.705, -3.979, -0.226,  # petal rotations in degrees
                    5.46, 167, 0, 12.51,  # petal translation and mag
                    -173.9, 36.1, -33.4, -43.6, -28.6, 152.3])  # FIF rotations
    bounds1 = ((-180, 180), (-180, 180), (-180, 180),  # in degrees
               (-300, 300), (-300, 300), (-1, 1), (12.4, 12.6),
               (-180, 180), (-180, 180), (-180, 180),
               (-180, 180), (-180, 180), (-180, 180))
    # assuming FIFs are perfectly straight, no bending
    # not allow in-plane translation here because it adds 12 free params
    solution1 = minimize(
        residue_ptl_fif, iv1, args=fif_data, bounds=bounds1, method='SLSQP',
        options={'disp': True, 'maxiter': 1e3, 'ftol': 1e-8})
    # what actual in-plane coherent shifts can cause the observed patterns
    # with all above fitted parameters fixed?
    iv2 = np.zeros(2 * len(fif_info))
    bounds2 = tuple([(-0.5, 0.5)] * len(iv2))
    solution2 = minimize(  # now allowing pinholes in-plane coherent shift
        residue_inplane_shifts, iv2, args=(solution1.x, fif_data),
        bounds=bounds2, method='SLSQP',
        options={'disp': True, 'maxiter': 1e3, 'ftol': 1e-9})
    modelprime, fif_ptl = plot_solutions(solution1, solution2)
    # path = os.path.join(data_dir, 'petal3.json')  # read existing json data
    # with open(path) as handle:
    #     data = json.load(handle)
    dump = {}  # write synthetic metrology as dict to json
    errors = (np.linalg.norm(fif_data - modelprime[:2, :], axis=0)
              / solution1.x[6])  # devide by magnification (conversion factor)
    for i, (fifid, info) in enumerate(fif_info.items()):
        dump[fifid] = {'petal_id': 0,
                       'device_loc': info[0]}
        for j, key in enumerate([
                'pinhole1', 'pinhole2', 'pinhole3', 'pinhole4', 'center']):
            dump[fifid][key] = {'x': fif_ptl[fifid][0, j],
                                'y': fif_ptl[fifid][1, j],
                                'z': fif_ptl[fifid][2, j]}
            if key == 'center':
                dump[fifid][key]['error'] = np.mean(errors[4*i:4*i+4])
            else:
                dump[fifid][key]['error'] = errors[4*i + j]
    with open(os.path.join(data_dir, 'petal0.json'), 'w') as handle:
        json.dump(dump, handle, ensure_ascii=False, sort_keys=False, indent=4)
