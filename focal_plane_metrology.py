# -*- coding: utf-8 -*-
"""
Created on Sun Jun  4 13:33:19 2017

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
from multiprocessing import Pool
import numpy as np
import pandas as pd
from scipy.optimize import minimize, minimize_scalar
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages

# 2017-11-22 (run 4)
petal_locations = [0, 1, 2, 3, 4, 6, 7, 8]  # lo of petals installed
petal_id_lookup = {0: '06',  # map between petal location and petal ID
                   1: '03',
                   2: '00',
                   3: '04',
                   4: '02',
                   5: '10',
                   6: '05',
                   7: '01',
                   8: '07',
                   9: '09'}

# # final alignment, 10 official petals
# petal_locations = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]  # lo of petals installed
# petal_id_lookup = {0: '04',  # map between petal location and petal ID
#                    1: '05',
#                    2: '06',
#                    3: '07',
#                    4: '08',
#                    5: '10',
#                    6: '11',
#                    7: '02',
#                    8: '03',
#                    9: '09'}

petal_ids = ['01', '02', '04', '00', '03', '05', '06', '07', '08', '09', '10',
             '11']  # petal production sequential order
fig_save_dir = r'C:\Users\givoltage\Google Drive\DESI\model_drawings\DESI Focal Plate Assy and Integration\FP Structure\metrology\duan_plots_and_data'
fig_save_dir = r'D:\20171122 (run 4)'
''' CMM data

2017/06
# 14 ft lbs run

4.4781,	-25.60003,	-83.09296
-107.04049,	-412.34349,	-102.98927
40.8575,	-423.96552,	-102.94535

# run 1
# PTL 01
4.4773,	-25.5948,	-83.12058
-107.04192,	-412.34,	-102.99234
40.85639,	-423.9622,	-102.94688

# run 2
# PTL01
4.45612,	-25.59908,	-83.11410
-107.38446,	-412.25313,	-102.98851
40.50410,	-423.99704,	-102.94503

# run 3
# PTL01
4.46046,	-25.59862,	-83.10798
-107.33292,	-412.26666,	-102.98560
40.55794,	-423.99204,	-102.94328

# run 4
# PTL01
4.45294,	-25.47059,	-84.71993
-107.63671,	-412.12974,	-103.06178
40.23908,	-423.96850,	-102.99510

# run 5
# PTL01
4.40935,	-25.55869,	-83.35969
-107.69158,	-412.14867,	-103.02441
40.19443,	-423.99777,	-102.96807

# PTL02
22.82853,	12.03275,	-83.54088
425.04125,	25.66662,	-103.54869
390.20098,	170.08625,	-102.84353

# PTL04
-25.71718,	3.66485,	-83.57341
-359.01328,	229.15198,	-103.09733
-415.87795,	91.82995,	-103.02373

2017/08/09 Petal1_2_3_4_5 Position0_4_6_8_2
dtb_pos_cmm = {0: np.array([[4.4270, -25.6063,  -83.2393],   # DTB 0
                            [-107.2668,  -412.3038, -102.9977],   # DTB 1
                            [40.6288,  -423.9956, -102.9422]]), # DTB 2
               1: np.array([[ 23.13790,  12.17047,  -82.36111],   # DTB 0
                            [425.43710,  24.27891, -102.92469],   # DTB 1
                            [391.25734, 168.64819, -102.90111]]), # DTB 2
               2: np.array([[25.7355,  -3.9524, -81.7074],   # DTB 0
                            [358.5042, -229.9112, -102.9793],   # DTB 1
                            [415.8817, -92.6059,  -102.9030]]), # DTB 2
               3: np.array([[22.82853,	12.03275,	-83.54088],   # DTB 0
                            [425.04125,	25.66662,	-103.54869],   # DTB 1
                            [390.20098,	170.08625,	-102.84353]]), # DTB 2
               4: np.array([[11.5034,  23.2342, -82.0262],   # DTB 0
                            [329.2688, 270.0592,  -103.4826],   # DTB 1
                            [216.3584, 366.6220,  -102.7970]]), # DTB 2
               5: np.array([[ 23.13790,  12.17047,  -82.36447],   # DTB 0
                            [425.43710,  24.27891, -102.91776],   # DTB 1
                            [391.25734, 168.64819, -102.89321]]), # DTB 2
               6: np.array([[-18.6401, 18.2689, -81.4708],   # DTB 0
                            [-155.0093,  396.7549,  -102.9882],   # DTB 1
                            [-281.8764,  319.2967 , -102.9671]]), # DTB 2
               7: np.array([[-25.71718,	3.66485,	-83.57341],   # DTB 0
                            [-359.01328,	229.15198,	-103.09733],   # DTB 1
                            [-415.87795,	91.82995,	-103.02373]]), # DTB 2
               8: np.array([[-23.1968, -12.1067,  -82.2729],   # DTB 0
                            [-425.0534,  -25.3937,  -102.9956],   # DTB 1
                            [-390.5567,  -169.9256, -103.0188]]), # DTB 2
               9: np.array([[ 23.13790,  12.17047,  -82.36447],   # DTB 0
                            [425.43710,  24.27891, -102.91776],   # DTB 1
                            [391.25734, 168.64819, -102.89321]]), # DTB 2
              }

2017/08/09 Petal1_2_3_4_5 Position0_4_6_8_2_A
dtb_pos_cmm = {0: np.array([[4.1872,	-25.7208,	-81.6723],   # DTB 0
                            [-107.2467,	-412.2113,	-102.9762],   # DTB 1
                            [41.0895,	-424.1034,	-102.8918]]), # DTB 2
               1: np.array([[ 23.13790,  12.17047,  -82.36111],   # DTB 0
                            [425.43710,  24.27891, -102.92469],   # DTB 1
                            [391.25734, 168.64819, -102.90111]]), # DTB 2
               2: np.array([[25.8319,	-3.8208,	-82.3193],   # DTB 0
                            [358.8527,	-229.1334,	-102.9946],   # DTB 1
                            [415.8386,	-91.9034,	-103.0219]]), # DTB 2
               3: np.array([[22.82853,	12.03275,	-83.54088],   # DTB 0
                            [425.04125,	25.66662,	-103.54869],   # DTB 1
                            [390.20098,	170.08625,	-102.84353]]), # DTB 2
               4: np.array([[11.4151,	23.3294,	-83.1942],   # DTB 0
                            [328.8552,	270.7950,	-103.0001],   # DTB 1
                            [215.9952,	367.0882,	-102.9447]]), # DTB 2
               5: np.array([[ 23.13790,  12.17047,  -82.36447],   # DTB 0
                            [425.43710,  24.27891, -102.91776],   # DTB 1
                            [391.25734, 168.64819, -102.89321]]), # DTB 2
               6: np.array([[-18.5823,	18.0751,	-82.0959],   # DTB 0
                            [-155.2874,	396.5124,	-103.4683],   # DTB 1
                            [-281.9814,	318.9159,	-102.7907]]), # DTB 2
               7: np.array([[-25.71718,	3.66485,	-83.57341],   # DTB 0
                            [-359.01328,	229.15198,	-103.09733],   # DTB 1
                            [-415.87795,	91.82995,	-103.02373]]), # DTB 2
               8: np.array([[-23.1689,	-12.1347,	-81.4020],   # DTB 0
                            [-425.2396,	-25.6328,	-102.9943],   # DTB 1
                            [-390.5036,	-170.1635,	-102.9757]]), # DTB 2
               9: np.array([[ 23.13790,  12.17047,  -82.36447],   # DTB 0
                            [425.43710,  24.27891, -102.91776],   # DTB 1
                            [391.25734, 168.64819, -102.89321]]), # DTB 2
              }

2017/08/09 Petal1_2_3_4_5 Position0_4_6_8_2_A_Center_Ring
dtb_pos_cmm = {0: np.array([[4.1921,	-25.7109,	-81.7634],   # DTB 0
                            [-107.2516,	-412.2183,	-103.0022],   # DTB 1
                            [41.0836,	-424.1112,	-102.9165]]), # DTB 2
               1: np.array([[ 23.13790,  12.17047,  -82.36111],   # DTB 0
                            [425.43710,  24.27891, -102.92469],   # DTB 1
                            [391.25734, 168.64819, -102.90111]]), # DTB 2
               2: np.array([[25.8265,	-3.8083,	-82.3541],   # DTB 0
                            [358.8424,	-229.1402,	-102.9886],   # DTB 1
                            [415.8300,	-91.9114,	-103.0040]]), # DTB 2
               3: np.array([[22.82853,	12.03275,	-83.54088],   # DTB 0
                            [425.04125,	25.66662,	-103.54869],   # DTB 1
                            [390.20098,	170.08625,	-102.84353]]), # DTB 2
               4: np.array([[11.4455,	23.3725,	-82.8678],   # DTB 0
                            [328.8623,	270.7992,	-102.9838],   # DTB 1
                            [216.0007,	367.0942,	-102.9391]]), # DTB 2
               5: np.array([[ 23.13790,  12.17047,  -82.36447],   # DTB 0
                            [425.43710,  24.27891, -102.91776],   # DTB 1
                            [391.25734, 168.64819, -102.89321]]), # DTB 2
               6: np.array([[-18.5789,	18.1010,	-81.9744],   # DTB 0
                            [-155.2896,	396.5123,	-103.4676],   # DTB 1
                            [-281.9832,	318.9133,	-102.7874]]), # DTB 2
               7: np.array([[-25.71718,	3.66485,	-83.57341],   # DTB 0
                            [-359.01328,	229.15198,	-103.09733],   # DTB 1
                            [-415.87795,	91.82995,	-103.02373]]), # DTB 2
               8: np.array([[-23.1181,	-12.1256,	-81.7091],   # DTB 0
                            [-425.2315,	-25.6386,	-102.9947],   # DTB 1
                            [-390.4957,	-170.1678,	-102.9800]]), # DTB 2
               9: np.array([[ 23.13790,  12.17047,  -82.36447],   # DTB 0
                            [425.43710,  24.27891, -102.91776],   # DTB 1
                            [391.25734, 168.64819, -102.89321]]), # DTB 2
              }

2017/08/09 Petal1_2_3_4_5 Position0_4_6_8_2_Center_Ring
dtb_pos_cmm = {0: np.array([[4.4375,	-25.6494,	-82.9543],   # DTB 0
                            [-107.2724,	-412.3054,	-103.0165],   # DTB 1
                            [40.6252,	-424.0028,	-102.9418]]), # DTB 2
               1: np.array([[ 23.13790,  12.17047,  -82.36111],   # DTB 0
                            [425.43710,  24.27891, -102.92469],   # DTB 1
                            [391.25734, 168.64819, -102.90111]]), # DTB 2
               2: np.array([[25.6946,	-3.9457,	-82.0011],   # DTB 0
                            [358.4965,	-229.9100,	-102.9356],   # DTB 1
                            [415.8752,	-92.6053,	-102.8695]]), # DTB 2
               3: np.array([[22.82853,	12.03275,	-83.54088],   # DTB 0
                            [425.04125,	25.66662,	-103.54869],   # DTB 1
                            [390.20098,	170.08625,	-102.84353]]), # DTB 2
               4: np.array([[11.5420,	23.2737,	-81.7051],   # DTB 0
                            [329.2849,	270.0662,	-103.5248],   # DTB 1
                            [216.3745,	366.6306,	-102.8486]]), # DTB 2
               5: np.array([[ 23.13790,  12.17047,  -82.36447],   # DTB 0
                            [425.43710,  24.27891, -102.91776],   # DTB 1
                            [391.25734, 168.64819, -102.89321]]), # DTB 2
               6: np.array([[-18.6197,	18.2013,	-81.9046],   # DTB 0
                            [-154.9941,	396.7444,	-102.9752],   # DTB 1
                            [-281.8636,	319.2927,	-102.9304]]), # DTB 2
               7: np.array([[-25.71718,	3.66485,	-83.57341],   # DTB 0
                            [-359.01328,	229.15198,	-103.09733],   # DTB 1
                            [-415.87795,	91.82995,	-103.02373]]), # DTB 2
               8: np.array([[-23.2258,	-12.1098,	-82.1010],   # DTB 0
                            [-425.0643,	-25.3830,	-102.9931],   # DTB 1
                            [-390.5664,	-169.9151,	-103.0403]]), # DTB 2
               9: np.array([[ 23.13790,  12.17047,  -82.36447],   # DTB 0
                            [425.43710,  24.27891, -102.91776],   # DTB 1
                            [391.25734, 168.64819, -102.89321]]), # DTB 2
              }

2017/08/09 Petal1_2_3_4_5 Position0_6_4_2_8_B
dtb_pos_cmm = {0: np.array([[4.3125,	-25.7488,	-81.4436],   # DTB 0
                            [-107.6097,	-412.1692,	-102.9695],   # DTB 1
                            [40.5679,	-423.9876,	-102.9520]]), # DTB 2
               1: np.array([[ 23.13790,  12.17047,  -82.36111],   # DTB 0
                            [425.43710,  24.27891, -102.92469],   # DTB 1
                            [391.25734, 168.64819, -102.90111]]), # DTB 2
               2: np.array([[25.6611,	-3.7420,	-82.0767],   # DTB 0
                            [358.5348,	-229.7812,	-103.4887],   # DTB 1
                            [415.5226,	-92.5762,	-102.8118]]), # DTB 2
               3: np.array([[22.82853,	12.03275,	-83.54088],   # DTB 0
                            [425.04125,	25.66662,	-103.54869],   # DTB 1
                            [390.20098,	170.08625,	-102.84353]]), # DTB 2
               4: np.array([[11.7554,	23.2369,	-81.6233],   # DTB 0
                            [329.4846,	269.8851,	-102.9853],   # DTB 1
                            [216.6307,	366.8838,	-102.8956]]), # DTB 2
               5: np.array([[ 23.13790,  12.17047,  -82.36447],   # DTB 0
                            [425.43710,  24.27891, -102.91776],   # DTB 1
                            [391.25734, 168.64819, -102.89321]]), # DTB 2
               6: np.array([[-18.6733,	18.2631,	-82.3087],   # DTB 0
                            [-155.5033,	396.3449,	-102.9803],   # DTB 1
                            [-282.2959,	318.8688,	-103.0141]]), # DTB 2
               7: np.array([[-25.71718,	3.66485,	-83.57341],   # DTB 0
                            [-359.01328,	229.15198,	-103.09733],   # DTB 1
                            [-415.87795,	91.82995,	-103.02373]]), # DTB 2
               8: np.array([[-23.0191,	-12.1486,	-83.1843],   # DTB 0
                            [-425.3163,	-24.9167,	-103.0161],   # DTB 1
                            [-390.9161,	-169.2275,	-102.9588]]), # DTB 2
               9: np.array([[ 23.13790,  12.17047,  -82.36447],   # DTB 0
                            [425.43710,  24.27891, -102.91776],   # DTB 1
                            [391.25734, 168.64819, -102.89321]]), # DTB 2
              }

2017/08/09 Petal1_2_3_4_5 Position0_6_4_2_8_B_Center_Ring

0 01  4.306151714 -25.65967876  -81.99537817
    -107.5974998  -412.1497337  -102.9052737
    40.58113183 -423.9656425  -102.8851082
2 00  25.71742726 -3.752707116  -81.75653534
    358.5624951 -229.7844672  -103.5367403
    415.5495519 -92.57838773  -102.8862439
4 02  11.70441665 23.18940413 -82.03236575
    329.4829831 269.8730733 -102.9891035
    216.6299119 366.8727973 -102.857263
6 04  -18.68022968  18.27740866 -82.20842223
    -155.5026633  396.3509597 -102.9421877
    -282.2944514  318.8704821 -103.0078759
8 03  -23.10723901  -12.18522014  -82.58056917
    -425.3336435  -24.9240861 -103.0534955
    -390.9336401  -169.2380798  -102.9815755

2017/08/09 Petal1_2_3_4_5 Position0_2_4_6_8_C_Center_Ring
dtb_pos_cmm = {0: np.array([[4.4409, -25.6956,  -82.1555],   # DTB 0
                            [-106.9750,  -412.4034, -102.9504],   # DTB 1
                            [40.9380,  -423.9888, -102.9159]]), # DTB 2
               1: np.array([[ 23.13790,  12.17047,  -82.36111],   # DTB 0
                            [425.43710,  24.27891, -102.92469],   # DTB 1
                            [391.25734, 168.64819, -102.90111]]), # DTB 2
               2: np.array([[25.7327,  -3.9070, -81.8177],   # DTB 0
                            [358.9055, -229.2893, -102.9819],   # DTB 1
                            [416.0439, -91.8863,  -102.9071]]), # DTB 2
               3: np.array([[22.82853,	12.03275,	-83.54088],   # DTB 0
                            [425.04125,	25.66662,	-103.54869],   # DTB 1
                            [390.20098,	170.08625,	-102.84353]]), # DTB 2
               4: np.array([[11.4818,  23.2879, -81.5309],   # DTB 0
                            [328.7666, 270.6905,  -103.4652],   # DTB 1
                            [215.6672, 367.0395,  -102.7877]]), # DTB 2
               5: np.array([[ 23.13790,  12.17047,  -82.36447],   # DTB 0
                            [425.43710,  24.27891, -102.91776],   # DTB 1
                            [391.25734, 168.64819, -102.89321]]), # DTB 2
               6: np.array([[-18.6736, 18.2012, -81.7480],   # DTB 0
                            [-155.9486,  396.3813,  -102.9985],   # DTB 1
                            [-282.6288,  318.6247 , -102.9751]]), # DTB 2
               7: np.array([[-25.71718,	3.66485,	-83.57341],   # DTB 0
                            [-359.01328,	229.15198,	-103.09733],   # DTB 1
                            [-415.87795,	91.82995,	-103.02373]]), # DTB 2
               8: np.array([[-23.2176, -12.1309,  -81.8876],   # DTB 0
                            [-425.0558,  -25.6160,  -102.9812],   # DTB 1
                            [-390.4803,  -170.1315, -103.0081]]), # DTB 2
               9: np.array([[ 23.13790,  12.17047,  -82.36447],   # DTB 0
                            [425.43710,  24.27891, -102.91776],   # DTB 1
                            [391.25734, 168.64819, -102.89321]]), # DTB 2
              }

# 2017-11-17 run 1
dtb_pos_cmm = {0: np.array([[4.1735,	-25.6740,	-81.9680],   # DTB 0
                            [-107.4897,	-412.1371,	-102.8910],   # DTB 1
                            [40.8334,	-424.1103,	-102.8359]]), # DTB 2
               1: np.array([[18.6743,	-18.1074,	-82.6164],   # DTB 0
                            [155.4888,	-396.6143,	-102.8889],   # DTB 1
                            [282.0296,	-319.1703,	-102.8710]]), # DTB 2
               2: np.array([[25.7062,	-3.6743,	-81.7622],   # DTB 0
                            [358.7403,	-229.4145,	-103.4131],   # DTB 1
                            [415.5979,	-92.1552,	-102.7286]]), # DTB 2
               3: np.array([[23.1556,	12.1231,	-82.0438],   # DTB 0
                            [425.0321,	25.1929,	-102.9075],   # DTB 1
                            [390.6092,	169.7450,	-102.9146]]), # DTB 2
               4: np.array([[11.7554,	23.2369,	-81.6233],   # DTB 0
                            [329.4846,	269.8851,	-102.9853],   # DTB 1
                            [216.6307,	366.8838,	-102.8956]]), # DTB 2
               5: np.array([[-4.3501,	25.7340,	-81.9638],   # DTB 0
                            [107.3870,	412.2233,	-102.9615],   # DTB 1
                            [-40.7836,	423.9673,	-102.9348]]), # DTB 2
               6: np.array([[-18.4454,	18.2190,	-82.0118],   # DTB 0
                            [-155.5146,	396.8222,	-102.9193],   # DTB 1
                            [-282.2076,	319.1092,	-102.9615]]), # DTB 2
               7: np.array([[-25.9144,	3.6912,	-82.0610],   # DTB 0
                            [-358.9540,	229.6060,	-102.9805],   # DTB 1
                            [-416.0889,	92.0545,	-102.9369]]), # DTB 2
               8: np.array([[-23.3388,	-12.2638,	-81.9155],   # DTB 0
                            [-425.3304,	-25.3912,	-102.9180],   # DTB 1
                            [-390.8561,	-169.8411,	-102.9696]]), # DTB 2
               9: np.array([[ 23.13790,  12.17047,  -82.36447],   # DTB 0
                            [425.43710,  24.27891, -102.91776],   # DTB 1
                            [391.25734, 168.64819, -102.89321]]), # DTB 2
              }

# 20171122 (run 4)
dtb_pos_cmm = {0: np.array([[4.418, -25.781, -81.858],   # DTB 0
                            [-107.553, -412.325, -102.839],   # DTB 1
                            [40.924, -424.167, -102.811]]), # DTB 2
               1: np.array([[18.663, -18.308, -81.955],   # DTB 0
                            [155.414, -396.399, -102.924],   # DTB 1
                            [282.221, -318.943, -102.959]]), # DTB 2
               2: np.array([[25.761, -3.876, -81.836],   # DTB 0
                            [358.742, -229.547, -102.913],   # DTB 1
                            [415.987, -92.192, -102.837]]), # DTB 2
               3: np.array([[23.109, 12.129, -81.804],   # DTB 0
                            [425.231, 25.245, -102.929],   # DTB 1
                            [390.630, 169.802, -102.902]]), # DTB 2
               4: np.array([[11.487, 23.292, -81.590],   # DTB 0
                            [329.002, 270.423, -103.408],   # DTB 1
                            [215.990, 366.870, -102.735]]), # DTB 2
               5: np.array([[-4.3501,	25.7340,	-81.9638],   # DTB 0
                            [107.3870,	412.2233,	-102.9615],   # DTB 1
                            [-40.7836,	423.9673,	-102.9348]]), # DTB 2
               6: np.array([[-18.476, 18.237, -81.790],   # DTB 0
                            [-155.512, 396.829, -102.775],   # DTB 1
                            [-282.212, 319.127, -102.832]]), # DTB 2
               7: np.array([[-25.826, 3.645, -82.360],  # DTB 0
                            [-358.962, 229.482, -102.939], # DTB 1
                            [-415.831, 92.454, -102.910]]), # DTB 2
               8: np.array([[-23.354, -12.225, -81.847],   # DTB 0
                            [-425.341, -25.373, -102.810],   # DTB 1
                            [-390.863, -169.819, -102.860]]), # DTB 2
               9: np.array([[ 23.13790,  12.17047,  -82.36447],   # DTB 0
                            [425.43710,  24.27891, -102.91776],   # DTB 1
                            [391.25734, 168.64819, -102.89321]]), # DTB 2
              }

# 20171122 (run 3)
dtb_pos_cmm = {0: np.array([[4.419, -25.789, -81.842],   # DTB 0
                            [-107.554, -412.326, -102.845],   # DTB 1
                            [40.924, -424.167, -102.816]]), # DTB 2
               1: np.array([[18.651, -18.301, -82.032],   # DTB 0
                            [155.410, -396.398, -102.929],   # DTB 1
                            [282.218, -318.941, -102.962]]), # DTB 2
               2: np.array([[25.749, -3.880, -81.876],   # DTB 0
                            [358.738, -229.545, -102.914],   # DTB 1
                            [415.983, -92.191, -102.839]]), # DTB 2
               3: np.array([[23.121, 12.133, -81.721],   # DTB 0
                            [425.230, 25.245, -102.932],   # DTB 1
                            [390.630, 169.802, -102.907]]), # DTB 2
               4: np.array([[11.487, 23.302, -81.542],   # DTB 0
                            [329.001, 270.425, -103.412],   # DTB 1
                            [215.991, 366.873, -102.738]]), # DTB 2
               5: np.array([[-4.3501,	25.7340,	-81.9638],   # DTB 0
                            [107.3870,	412.2233,	-102.9615],   # DTB 1
                            [-40.7836,	423.9673,	-102.9348]]), # DTB 2
               6: np.array([[-18.470, 18.237, -81.804],   # DTB 0
                            [-155.511, 396.830, -102.774],   # DTB 1
                            [-282.211, 319.128, -102.834]]), # DTB 2
               7: np.array([[-25.825, 3.646, -82.346],   # DTB 0
                            [-358.961, 229.482, -102.942],   # DTB 1
                            [-415.829, 92.456, -102.914]]), # DTB 2
               8: np.array([[-23.357, -12.230, -81.820],   # DTB 0
                            [-425.340, -25.373, -102.816],   # DTB 1
                            [-390.862, -169.818, -102.865]]), # DTB 2
               9: np.array([[ 23.13790,  12.17047,  -82.36447],   # DTB 0
                            [425.43710,  24.27891, -102.91776],   # DTB 1
                            [391.25734, 168.64819, -102.89321]]), # DTB 2
              }

# 20171122 (run 2)
dtb_pos_cmm = {0: np.array([[4.362, -25.710, -81.810],   # DTB 0
                            [-107.290, -412.240, -102.917],   # DTB 1
                            [40.886, -423.953, -102.895]]), # DTB 2
               1: np.array([[18.450, -18.236, -81.783],   # DTB 0
                            [155.607, -396.788, -102.776],   # DTB 1
                            [282.282, -319.046, -102.833]]), # DTB 2
               2: np.array([[25.917, -3.711, -81.827],   # DTB 0
                            [358.954, -229.606, -102.861],   # DTB 1
                            [416.083, -92.050, -102.831]]), # DTB 2
               3: np.array([[23.328, 12.232, -81.836],   # DTB 0
                            [425.314, 25.397, -102.813],   # DTB 1
                            [390.827, 169.845, -102.865]]), # DTB 2
               4: np.array([[11.487, 23.302, -81.542],   # DTB 0
                            [329.001, 270.425, -103.412],   # DTB 1
                            [215.991, 366.873, -102.738]]), # DTB 2
               5: np.array([[-4.225, 25.692, -81.848],   # DTB 0
                            [107.444, 412.134, -102.905],   # DTB 1
                            [-40.878, 424.114, -102.825]]), # DTB 2
               6: np.array([[-18.734, 18.110, -82.381],   # DTB 0
                            [-155.544, 396.620, -102.928],   # DTB 1
                            [-282.086, 319.174, -102.903]]), # DTB 2
               7: np.array([[-25.724, 3.686, -81.623],   # DTB 0
                            [-358.778, 229.445, -103.418],   # DTB 1
                            [-415.642, 92.188, -102.750]]), # DTB 2
               8: np.array([[-23.228, -12.090, -81.956],   # DTB 0
                            [-425.072, -25.290, -102.936],   # DTB 1
                            [-390.603, -169.830, -102.961]]), # DTB 2
               9: np.array([[ 23.13790,  12.17047,  -82.36447],   # DTB 0
                            [425.43710,  24.27891, -102.91776],   # DTB 1
                            [391.25734, 168.64819, -102.89321]]), # DTB 2
              }

# 20171122 (run 1)
dtb_pos_cmm = {0: np.array([[4.363, -25.716, -81.781],   # DTB 0
                            [-107.290, -412.239, -102.913],   # DTB 1
                            [40.886, -423.952, -102.894]]), # DTB 2
               1: np.array([[18.453, -18.235, -81.778],   # DTB 0
                            [155.606, -396.784, -102.778],   # DTB 1
                            [282.281, -319.044, -102.837]]), # DTB 2
               2: np.array([[25.910, -3.708, -81.857],   # DTB 0
                            [358.952, -229.604, -102.866],   # DTB 1
                            [416.081, -92.049, -102.834]]), # DTB 2
               3: np.array([[23.333, 12.232, -81.811],   # DTB 0
                            [425.312, 25.398, -102.815],   # DTB 1
                            [390.825, 169.845, -102.864]]), # DTB 2
               4: np.array([[11.487, 23.302, -81.542],   # DTB 0
                            [329.001, 270.425, -103.412],   # DTB 1
                            [215.991, 366.873, -102.738]]), # DTB 2
               5: np.array([[-4.228, 25.694, -81.844],   # DTB 0
                            [107.443, 412.133, -102.903],   # DTB 1
                            [-40.879, 424.113, -102.826]]), # DTB 2
               6: np.array([[-18.736, 18.116, -82.352],   # DTB 0
                            [-155.545, 396.622, -102.931],   # DTB 1
                            [-282.086, 319.175, -102.906]]), # DTB 2
               7: np.array([[-25.731, 3.691, -81.579],   # DTB 0
                            [-358.778, 229.447, -103.419],   # DTB 1
                            [-415.642, 92.188, -102.749]]), # DTB 2
               8: np.array([[-23.223, -12.084, -81.984],   # DTB 0
                            [-425.070, -25.289, -102.934],   # DTB 1
                            [-390.601, -169.828, -102.956]]), # DTB 2
               9: np.array([[ 23.13790,  12.17047,  -82.36447],   # DTB 0
                            [425.43710,  24.27891, -102.91776],   # DTB 1
                            [391.25734, 168.64819, -102.89321]]), # DTB 2
              }

"2018-01-25 (run 4)"
dtb_pos_cmm = {0: np.array([[4.362, -25.713, -81.798],   # DTB 0
                            [-107.414, -412.211, -102.913],   # DTB 1
                            [40.760, -423.970, -102.891]]), # DTB 2
               1: np.array([[18.447, -18.242, -81.779],   # DTB 0
                            [155.481, -396.838, -102.773],   # DTB 1
                            [282.181, -319.136, -102.831]]), # DTB 2
               2: np.array([[25.917, -3.717, -81.825],   # DTB 0
                            [358.885, -229.714, -102.860],   # DTB 1
                            [416.055, -92.175, -102.830]]), # DTB 2
               3: np.array([[23.338, 12.239, -81.776],   # DTB 0
                            [425.317, 25.412, -102.809],   # DTB 1
                            [390.827, 169.856, -102.861]]), # DTB 2
               4: np.array([[11.671, 23.542, -81.893],   # DTB 0
                            [329.177, 270.374, -102.877],   # DTB 1
                            [216.397, 367.140, -102.898]]), # DTB 2
               5: np.array([[-4.388, 25.952, -81.701],   # DTB 0
                            [107.598, 412.227, -102.772],   # DTB 1
                            [-40.727, 424.155, -102.756]]), # DTB 2
               6: np.array([[-18.723, 18.200, -81.796],   # DTB 0
                            [-155.290, 396.507, -102.851],   # DTB 1
                            [-282.285, 318.730, -102.672]]), # DTB 2
               7: np.array([[-25.726, 3.692, -81.606],   # DTB 0
                            [-358.859, 229.332, -103.417],   # DTB 1
                            [-415.675, 92.053, -102.750]]), # DTB 2
               8: np.array([[-23.230, -12.090, -81.938],   # DTB 0
                            [-425.072, -25.278, -102.934],   # DTB 1
                            [-390.606, -169.821, -102.957]]), # DTB 2
               9: np.array([[-11.496, -23.316,  -81.685],   # DTB 0
                            [-329.255, -270.513, -102.869],   # DTB 1
                            [-216.270, -367.187, -102.845]]), # DTB 2
              }

'''


# DTB coordinates in CS5 measured by CMM
# 2018-01-25 (run 4)
# 20171122 (run 4)
dtb_pos_cmm = {0: np.array([[4.418, -25.781, -81.858],   # DTB 0
                            [-107.553, -412.325, -102.839],   # DTB 1
                            [40.924, -424.167, -102.811]]), # DTB 2
               1: np.array([[18.663, -18.308, -81.955],   # DTB 0
                            [155.414, -396.399, -102.924],   # DTB 1
                            [282.221, -318.943, -102.959]]), # DTB 2
               2: np.array([[25.761, -3.876, -81.836],   # DTB 0
                            [358.742, -229.547, -102.913],   # DTB 1
                            [415.987, -92.192, -102.837]]), # DTB 2
               3: np.array([[23.109, 12.129, -81.804],   # DTB 0
                            [425.231, 25.245, -102.929],   # DTB 1
                            [390.630, 169.802, -102.902]]), # DTB 2
               4: np.array([[11.487, 23.292, -81.590],   # DTB 0
                            [329.002, 270.423, -103.408],   # DTB 1
                            [215.990, 366.870, -102.735]]), # DTB 2
               5: np.array([[-4.3501, 25.7340,  -81.9638],   # DTB 0
                            [107.3870,  412.2233, -102.9615],   # DTB 1
                            [-40.7836,  423.9673, -102.9348]]), # DTB 2
               6: np.array([[-18.476, 18.237, -81.790],   # DTB 0
                            [-155.512, 396.829, -102.775],   # DTB 1
                            [-282.212, 319.127, -102.832]]), # DTB 2
               7: np.array([[-25.826, 3.645, -82.360],  # DTB 0
                            [-358.962, 229.482, -102.939], # DTB 1
                            [-415.831, 92.454, -102.910]]), # DTB 2
               8: np.array([[-23.354, -12.225, -81.847],   # DTB 0
                            [-425.341, -25.373, -102.810],   # DTB 1
                            [-390.863, -169.819, -102.860]]), # DTB 2
               9: np.array([[ 23.13790,  12.17047,  -82.36447],   # DTB 0
                            [425.43710,  24.27891, -102.91776],   # DTB 1
                            [391.25734, 168.64819, -102.89321]]), # DTB 2
              }

# %% function definitions

'''
rotation matrices, input in radians

This is a rotation whose yaw, pitch, and roll angles are α, β and γ,
or more formally intrinsic rotation whose Tait-Bryan angles are α, β, γ
about axes z, y, x respectively.

For intrinsic rotations, the coordinate system is rotated.

'''


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


def R_general(axis, angle):
    '''
    axis must be a 3-element unit vector

    the rotation formalism is outlined in
    http://ksuweb.kennesaw.edu/~plaval/math4490/rotgen.pdf
    where the arbitrary axis passes through the origin.

    '''
    ux = axis[0]
    uy = axis[1]
    uz = axis[2]
    C = np.cos(angle)
    S = np.sin(angle)
    t = 1-C
    R_general = np.array([
                          [t*ux**2+C,    t*ux*uy-S*uz, t*ux*uz+S*uy],
                          [t*ux*uy+S*uz, t*uy**2+C,    t*uy*uz-S*ux],
                          [t*ux*uz-S*uy, t*uy*uz+S*ux, t*uz**2+C]])
    return R_general


def angles_to_unit_vector(theta_deg, phi_deg):
    # returns a unit vector from polar and azimuthal angles
    theta = np.radians(theta_deg)
    phi = np.radians(phi_deg)
    return np.array([np.sin(theta)*np.cos(phi),
                     np.sin(theta)*np.sin(phi),
                     np.cos(theta)])


def vector_to_angles(x, y, z):
    r = np.sqrt(np.square(x) + np.square(y) + np.square(z))
    theta = np.arccos(z/r)  # radians
    phi = np.arctan(y/x)
    return [np.degrees(theta), np.degrees(phi)]


def throughput_tilt(tilt):
    # takes degree input
    throughput_tilt = -0.0133*np.power(tilt, 2) - 0.0175*tilt + 1.0
    return np.multiply(throughput_tilt, throughput_tilt > 0)


def throughput_defocus(delta_f):
    delta_f_um = np.abs(delta_f) * 1000
    # takes micron input
    throughput_defocus = (- 1.804e-14*np.power(delta_f_um, 5)
                          + 1.593e-11*np.power(delta_f_um, 4)
                          - 5.955e-10*np.power(delta_f_um, 3)
                          - 3.433e-6*np.power(delta_f_um, 2)
                          + 3.251e-7*delta_f_um
                          + 1.0)
    return np.multiply(throughput_defocus, throughput_defocus > 0)


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


# %% evaluate petal throughput and calculate optimal 1DF alignment

def evaluate_petal(petal_location):

    petal_id = petal_id_lookup[petal_location]
    # read in dataframes
    df = pd.read_pickle(os.path.join(fig_save_dir,
                                     str(petal_ids.index(petal_id)+1).zfill(2)
                                     + '-ptl_' + petal_id + '-df_data.pickle'),
                        compression='gzip')

    df_transformations = pd.read_pickle(os.path.join(
            fig_save_dir,
            str(petal_ids.index(petal_id)+1).zfill(2)
            + '-ptl_' + petal_id + '-df_transformations.pickle'),
        compression='gzip')

    # %% from 3 DTB positions, calculate all petal info and throughput

    df.loc[('diameter', 'ACT')] = df.loc[('diameter', 'ABC')].values
    # create arrays for datum tooling ball positions in two alignments
    dtb_pos_abc = np.concatenate((
            df.loc[('dtb_x', 'ABC'), 'actual'].values.reshape(-1, 1),
            df.loc[('dtb_y', 'ABC'), 'actual'].values.reshape(-1, 1),
            df.loc[('dtb_z', 'ABC'), 'actual'].values.reshape(-1, 1)),
        axis=1).T.astype(np.float64)[:, :3]
    dtb_pos_act_cs5 = dtb_pos_cmm[petal_location].T
    dtb_pos_act = np.empty(dtb_pos_act_cs5.shape)
    # the measured actual positions are in CS5
    # need to rotate by an integer multiple of 36 degrees
    for j in range(3):
        dtb_pos_act[:, j] = cs5_to_petal(dtb_pos_act_cs5[:, j], petal_location)
        df.loc[('dtb_x', 'ACT'), 'actual'][j] = dtb_pos_act[0, j]
        df.loc[('dtb_y', 'ACT'), 'actual'][j] = dtb_pos_act[1, j]
        df.loc[('dtb_z', 'ACT'), 'actual'][j] = dtb_pos_act[2, j]

    def total_residue(parameters):
        # this is the function to be minimised
        alpha = parameters[0]
        beta = parameters[1]
        gamma = parameters[2]
        T = parameters[3:]
        R = Rxyz(alpha, beta, gamma)
        dtb_pos_abc_rot = np.empty(dtb_pos_abc.shape)  # rotated from ABC
        # rotate each column of x_abc and fill x_abc_rot
        for j in range(dtb_pos_abc.shape[1]):
            dtb_pos_abc_rot[:, j] = matmul(R, dtb_pos_abc[:, j]) + T

        return np.sum(np.linalg.norm(dtb_pos_abc_rot - dtb_pos_act), axis=0)

    # minimisation routine
    # p0 = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    p0 = df_transformations.loc['ZBF'].values
    bounds = ((-np.pi/2, np.pi//2), (-np.pi/2, np.pi/2), (-np.pi/2, np.pi/2),
              (-10, 10), (-10, 10), (-10, 10))
    solution_act = minimize(total_residue, p0,
                            bounds=bounds,
                            method='SLSQP',
                            options={'disp': True,
                                     'maxiter': 1000})
    print('PTL'+petal_id+' ACT transformation found as \n'
          + 'X Rotation (Roll) : {}° \n'.format(np.degrees(solution_act.x[0]))
          + 'Y Rotation (Pitch): {}° \n'.format(np.degrees(solution_act.x[1]))
          + 'Z Rotation (Yaw)  : {}° \n'.format(np.degrees(solution_act.x[2]))
          + 'Translation: {} \n'.format(solution_act.x[3:])
          + 'With least square: {} \n'.format(solution_act.fun))

    # save transformation parameters
    df_transformations.loc['ACT'] = solution_act.x

    # Calculate R and T
    parameters = solution_act.x
    alpha = parameters[0]
    beta = parameters[1]
    gamma = parameters[2]
    T = parameters[3:]
    R = Rxyz(alpha, beta, gamma)  # yaw-pitch-roll system

    # write results with these parameters to dataframe
    x_abc = np.concatenate((
            df.loc[('x', 'ABC'), 'actual'].values.reshape(-1, 1),
            df.loc[('y', 'ABC'), 'actual'].values.reshape(-1, 1),
            df.loc[('z', 'ABC'), 'actual'].values.reshape(-1, 1)),
        axis=1).T.astype(np.float64)
    theta_abc = np.radians(df.loc[('nutation', 'ABC'), 'actual']
                           .values.astype(np.float64))
    phi_abc = np.radians(df.loc[('precession', 'ABC'), 'actual']
                         .values.astype(np.float64))
    x0_abc = x_abc - np.array([np.sin(theta_abc)*np.cos(phi_abc),
                               np.sin(theta_abc)*np.sin(phi_abc),
                               np.cos(theta_abc)])  # another point along axis

    x_abc_rot = np.empty(x_abc.shape)  # rotated matrix from ABC
    x0_abc_rot = np.empty(x_abc.shape)

    # rotate each column of x_abc and fill x_abc_rot
    for j in range(x_abc.shape[1]):
        x_abc_rot[:, j] = matmul(R, x_abc[:, j]) + T
        x0_abc_rot[:, j] = matmul(R, x0_abc[:, j]) + T

    n_rot = x_abc_rot - x0_abc_rot  # rotated axis direction
    [theta_rot, phi_rot] = vector_to_angles(
            n_rot[0, :], n_rot[1, :], n_rot[2, :])

    df.loc[('x', 'ACT'), 'actual'] = x_abc_rot[0, :]
    df.loc[('y', 'ACT'), 'actual'] = x_abc_rot[1, :]
    df.loc[('z', 'ACT'), 'actual'] = x_abc_rot[2, :]
    df.loc[('nutation', 'ACT'), 'actual'] = theta_rot
    df.loc[('precession', 'ACT'), 'actual'] = phi_rot
    df.loc[('x', 'ACT'), 'deviation'] = df.loc[('x', 'ACT'), 'actual'].values \
        - df.loc[('x', 'ACT'), 'nominal'].values
    df.loc[('y', 'ACT'), 'deviation'] = df.loc[('y', 'ACT'), 'actual'].values \
        - df.loc[('y', 'ACT'), 'nominal'].values
    df.loc[('z', 'ACT'), 'deviation'] = df.loc[('z', 'ACT'), 'actual'].values \
        - df.loc[('z', 'ACT'), 'nominal'].values
    df.loc[('nutation', 'ACT'), 'deviation'] = \
        df.loc[('nutation', 'ACT'), 'actual'].values \
        - df.loc[('nutation', 'ACT'), 'nominal'].values
    df.loc[('precession', 'ACT'), 'deviation'] = \
        df.loc[('precession', 'ACT'), 'actual'].values \
        - df.loc[('precession', 'ACT'), 'nominal'].values

    # calculate r
    df.loc[('r', 'ACT'), 'actual'] = np.sqrt(
        np.square(df.loc[('x', 'ACT'), 'actual'].values.astype(np.float64))
        + np.square(df.loc[('y', 'ACT'), 'actual'].values.astype(np.float64)))

    # calculate tilt
    theta0 = df.loc['nutation', 'ACT']['nominal'].values.astype(np.float64)
    theta = df.loc['nutation', 'ACT']['actual'].values.astype(np.float64)
    phi0 = df.loc['precession', 'ACT']['nominal'].values.astype(np.float64)
    phi = df.loc['precession', 'ACT']['actual'].values.astype(np.float64)

    n0 = angles_to_unit_vector(theta0, phi0)
    n = angles_to_unit_vector(theta, phi)
    tilt = np.array(
        [np.degrees(np.arccos(np.dot(n0[:, i], n[:, i]))) for i in range(514)])
    df.loc[('tilt', 'ACT'), 'actual'] = tilt

    # calculate defocus
    delta_r = np.array([
                df.loc[('x', 'ACT'), 'deviation'].values.astype(np.float64),
                df.loc[('y', 'ACT'), 'deviation'].values.astype(np.float64),
                df.loc[('z', 'ACT'), 'deviation'].values.astype(np.float64)
                ])
    delta_f = np.array([np.dot(delta_r[:, i], n0[:, i]) for i in range(514)])
    df.loc[('defocus', 'ACT'), 'actual'] = delta_f

    # calculate combined throughput
    throughput = throughput_tilt(tilt) * throughput_defocus(delta_f)
    df.loc[('throughput', 'ACT'), 'actual'] = throughput
    df.loc[('throughput', 'ACT'), 'deviation'] = 1 - throughput

    # %% 1DF alignment, based on 3 DTb positions
    '''
    the predicted 1 degree of freedom that is actually adjustable
    is an axis of rotation passing through two points:
        [438.99959, 0.6, -108]
        [355.1584605, 258.0377258, -108]
    '''

    df.loc[('diameter', '1DF')] = df.loc[('diameter', 'ABC')].values

    axp1 = np.array([355.1584605, 258.0377258, -108])
    axp2 = np.array([438.99959,   0.6,         -108])
    axis = (axp1-axp2)/np.linalg.norm(axp1-axp2)

    x_act = np.concatenate((
            df.loc[('x', 'ACT'), 'actual'].values.reshape(-1, 1),
            df.loc[('y', 'ACT'), 'actual'].values.reshape(-1, 1),
            df.loc[('z', 'ACT'), 'actual'].values.reshape(-1, 1)),
        axis=1).T.astype(np.float64)
    theta_act = np.radians(df.loc[('nutation', 'ACT'), 'actual']
                           .values.astype(np.float64))
    phi_act = np.radians(df.loc[('precession', 'ACT'), 'actual']
                           .values.astype(np.float64))
    x0_act = x_act - np.array([np.sin(theta_act)*np.cos(phi_act),
                               np.sin(theta_act)*np.sin(phi_act),
                               np.cos(theta_act)])  # another point along axis

    def throughput_loss_min_1df(angle):
        # this is the function to be minimised
        R = R_general(axis, angle)
        x_act_rot = np.empty(x_act.shape)  # rotated matrix from ABC
        x0_act_rot = np.empty(x_act.shape)
        # rotate each column of x_abc and fill x_abc_rot
        for k in range(x_act.shape[1]):
            x_act_rot[:, k] = matmul(R, x_act[:, k] - axp1) + axp1
            x0_act_rot[:, k] = matmul(R, x0_act[:, k] - axp1) + axp1
        n_rot = x_act_rot - x0_act_rot  # not necessarily of unit length
        [theta_rot, phi_rot] = vector_to_angles(
                n_rot[0, :], n_rot[1, :], n_rot[2, :])
        delta_r = np.array(
                [x_act_rot[0, :] - df.loc[('x', 'ABC'), 'nominal'].values,
                 x_act_rot[1, :] - df.loc[('y', 'ABC'), 'nominal'].values,
                 x_act_rot[2, :] - df.loc[('z', 'ABC'), 'nominal'].values])
        # calculate tilt
        n = angles_to_unit_vector(theta_rot, phi_rot)  # ensure unit length
        tilt = np.array(
                [np.degrees(np.arccos(np.dot(n0[:, i], n[:, i])))
                 for i in range(514)])
        delta_f = np.array(
                [np.dot(delta_r[:, i], n0[:, i]) for i in range(514)])
        throughput = throughput_tilt(tilt) * throughput_defocus(delta_f)

        return np.mean(1-throughput)

    solution_1df = minimize_scalar(throughput_loss_min_1df,
                                   bounds=(-np.pi/2, np.pi/2),
                                   method='Brent')

    print('PTL'+petal_id+' 1DF transformation found as \n'
          + 'Rotation: {}° \n'.format(np.degrees(solution_1df.x))
          + 'With least mean throughput loss: {} \n'.format(solution_1df.fun)
          + 'Compared with current throughput loss: {} \n'.format(
                  np.mean(df.loc[('throughput', 'ACT'), 'deviation'].values)))

    # save transformation parameters
    df_transformations.loc['1DF', 'alpha'] = solution_1df.x

    R = R_general(axis, solution_1df.x)
    x_act_rot = np.empty(x_act.shape)  # rotated matrix from ABC
    x0_act_rot = np.empty(x_act.shape)
    # rotate each column of x_abc and fill x_abc_rot
    for k in range(x_act.shape[1]):
        x_act_rot[:, k] = matmul(R, x_act[:, k] - axp1) + axp1
        x0_act_rot[:, k] = matmul(R, x0_act[:, k] - axp1) + axp1
    n_rot = x_act_rot - x0_act_rot  # not necessarily of unit length
    [theta_rot, phi_rot] = vector_to_angles(
            n_rot[0, :], n_rot[1, :], n_rot[2, :])

    df.loc[('x', '1DF'), 'actual'] = x_act_rot[0, :]
    df.loc[('y', '1DF'), 'actual'] = x_act_rot[1, :]
    df.loc[('z', '1DF'), 'actual'] = x_act_rot[2, :]
    df.loc[('nutation', '1DF'), 'actual'] = theta_rot
    df.loc[('precession', '1DF'), 'actual'] = phi_rot
    df.loc[('x', '1DF'), 'nominal'] = df.loc[('x', 'ABC'), 'nominal'].values
    df.loc[('y', '1DF'), 'nominal'] = df.loc[('y', 'ABC'), 'nominal'].values
    df.loc[('z', '1DF'), 'nominal'] = df.loc[('z', 'ABC'), 'nominal'].values
    df.loc[('nutation', '1DF'), 'nominal'] = \
        df.loc[('nutation', 'ABC'), 'nominal'].values
    df.loc[('precession', '1DF'), 'nominal'] = \
        df.loc[('precession', 'ABC'), 'nominal'].values
    df.loc[('x', '1DF'), 'deviation'] = \
        df.loc[('x', '1DF'), 'actual'].values \
        - df.loc[('x', '1DF'), 'nominal'].values
    df.loc[('y', '1DF'), 'deviation'] = \
        df.loc[('y', '1DF'), 'actual'].values \
        - df.loc[('y', '1DF'), 'nominal'].values
    df.loc[('z', '1DF'), 'deviation'] = \
        df.loc[('z', '1DF'), 'actual'].values \
        - df.loc[('z', '1DF'), 'nominal'].values
    df.loc[('nutation', '1DF'), 'deviation'] = \
        df.loc[('nutation', '1DF'), 'actual'].values \
        - df.loc[('nutation', '1DF'), 'nominal'].values
    df.loc[('precession', '1DF'), 'deviation'] = \
        df.loc[('precession', '1DF'), 'actual'].values \
        - df.loc[('precession', '1DF'), 'nominal'].values

    # datum tooilng balls
    for j in range(3):
        x = np.array([df.loc[('dtb_x', 'ACT'), 'actual'][j],
                      df.loc[('dtb_y', 'ACT'), 'actual'][j],
                      df.loc[('dtb_z', 'ACT'), 'actual'][j]])
        xp = matmul(R, x - axp1) + axp1
        df.loc[('dtb_x', '1DF'), 'actual'][j] = xp[0]
        df.loc[('dtb_y', '1DF'), 'actual'][j] = xp[1]
        df.loc[('dtb_z', '1DF'), 'actual'][j] = xp[2]

    # calculate r
    df.loc[('r', '1DF'), 'actual'] = np.sqrt(
        np.square(df.loc[('x', '1DF'), 'actual'].values.astype(np.float64))
        + np.square(df.loc[('y', '1DF'), 'actual'].values.astype(np.float64)))

    # calculate tilt
    theta0 = df.loc['nutation', 'ABC']['nominal'].values.astype(np.float64)
    theta = df.loc['nutation', '1DF']['actual'].values.astype(np.float64)
    phi0 = df.loc['precession', 'ABC']['nominal'].values.astype(np.float64)
    phi = df.loc['precession', '1DF']['actual'].values.astype(np.float64)

    n0 = angles_to_unit_vector(theta0, phi0)
    n = angles_to_unit_vector(theta, phi)
    tilt = np.array(
            [np.degrees(np.arccos(np.dot(n0[:, i], n[:, i])))
             for i in range(514)])
    df.loc[('tilt', '1DF'), 'actual'] = tilt

    # calculate defocus
    delta_r = np.array([
                df.loc[('x', '1DF'), 'deviation'].values.astype(np.float64),
                df.loc[('y', '1DF'), 'deviation'].values.astype(np.float64),
                df.loc[('z', '1DF'), 'deviation'].values.astype(np.float64)
                ])
    delta_f = np.array([np.dot(delta_r[:, i], n0[:, i]) for i in range(514)])
    df.loc[('defocus', '1DF'), 'actual'] = delta_f

    # calculate combined throughput
    throughput = throughput_tilt(tilt) * throughput_defocus(delta_f)
    df.loc[('throughput', '1DF'), 'actual'] = throughput
    df.loc[('throughput', '1DF'), 'deviation'] = 1 - throughput

    # %% all plots
    figtitle_prefix = 'Petal ' + str(petal_id)
    figtitles = {'diameter': 'Cylinder Diameter Deviation',
                 'x': 'Spotface Centre X Deviation',
                 'y': 'Spotface Centre Y Deviation',
                 'z': 'Spotface Centre Z Deviation',
                 'nutation': 'Nutation Deviation',
                 'precession': 'Precession Deviation',
                 'tilt': 'Cylinder Axial Tilt',
                 'defocus': 'Spotface Centre Defocus',
                 'throughput': 'Throughput Percentage'}
    axtitles = {'diameter': r'$\delta D/\mathrm{mm}$',
                'x': r'$\delta x/\mathrm{mm}$',
                'y': r'$\delta y/\mathrm{mm}$',
                'z': r'$\delta z/\mathrm{mm}$',
                'nutation': r'$\delta \theta/\degree$',
                'precession': r'$\delta \varphi/\degree$',
                'tilt': r'$\delta/\degree$',
                'defocus': r'$\delta f/\mathrm{mm}$',
                'throughput': r'$\eta \times 100\%$'}
#    colour_range = {'diameter': [0.008, 0.018],
#                    'x': [-0.03, 0.03],
#                    'y': [-0.03, 0.03],
#                    'z': [-0.03, 0.03],
#                    'nutation': [-0.028, 0.028],
#                    'precession': [-0.028, 0.028],
#                    'tilt': [0, 0.05],
#                    'defocus': [-0.03, 0.03],
#                    'throughput': [0, 0.005]}
    tol_lower = {'diameter': df.loc['diameter', 'ABC']['lowertol'],
                 'x': df.loc['x', 'ABC']['lowertol'],
                 'y': df.loc['y', 'ABC']['lowertol'],
                 'z': df.loc['z', 'ABC']['lowertol'],
                 'nutation': df.loc['nutation', 'ABC']['lowertol'],
                 'precession': df.loc['precession', 'ABC']['lowertol'],
                 'tilt': df.loc['tilt', 'ABC']['lowertol'],
                 'defocus': df.loc['defocus', 'ABC']['lowertol'],
                 'throughput': df.loc['throughput', 'ABC']['lowertol']*100}
    tol_upper = {'diameter': df.loc['diameter', 'ABC']['uppertol'],
                 'x': df.loc['x', 'ABC']['uppertol'],
                 'y': df.loc['y', 'ABC']['uppertol'],
                 'z': df.loc['z', 'ABC']['uppertol'],
                 'nutation': df.loc['nutation', 'ABC']['uppertol'],
                 'precession': df.loc['precession', 'ABC']['uppertol'],
                 'tilt': df.loc['tilt', 'ABC']['uppertol'],
                 'defocus': df.loc['defocus', 'ABC']['uppertol'],
                 'throughput': df.loc['throughput', 'ABC']['uppertol']*100}
    units = {'diameter': ' mm',
             'x': ' mm',
             'y': ' mm',
             'z': ' mm',
             'nutation': r'$\degree$',
             'precession': r'$\degree$',
             'tilt': r'$\degree$',
             'defocus': ' mm',
             'throughput': '%'}

    # for alignment in ['ACT', '1DF']:
    for alignment in ['ACT']:

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
                   'throughput': df.loc['throughput', alignment]['actual']*100}

        for feature in ['diameter', 'x', 'y', 'z', 'nutation', 'precession',
                        'tilt', 'defocus', 'throughput']:

            fig = plt.figure(figsize=(18, 6))
            gs = gridspec.GridSpec(1, 2, width_ratios=[3, 2])
            fig.suptitle(figtitle_prefix + ' ' + figtitles[feature] + ' ('
                         + alignment + ' Alignment)')
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
            plot0 = ax0.scatter(x, y, marker='o', lw=4,
                                c=colours[feature], cmap='plasma')
            # vmin = colour_range[feature][0], vmax = colour_range[feature][1]
            fig.colorbar(plot0, ax=ax0)
            ax0.axis('equal')
            ax0.set_title(axtitles[feature])
            ax0.set_xlabel('x/mm')
            ax0.set_ylabel('y/mm')
            ax0.text(50, 150, textstr, fontsize=12, linespacing=1.5,
                     bbox=textbbox, ha='left', va='bottom')
            ax1 = plt.subplot(gs[1])
            ax1.plot(r, colours[feature], 'b.', r, tol_upper[feature], 'r--',
                     r, tol_lower[feature], 'r--', )
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

    # export dataframes
    df.to_csv(os.path.join(fig_save_dir,
                           str(petal_ids.index(petal_id)+1).zfill(2)
                           + '-ptl_' + petal_id + '-df_data.csv'))
    df.to_pickle(os.path.join(fig_save_dir,
                              str(petal_ids.index(petal_id)+1).zfill(2)
                              + '-ptl_' + petal_id + '-df_data.pickle'),
                 compression='gzip')
    df_transformations.to_csv(os.path.join(
            fig_save_dir,
            str(petal_ids.index(petal_id)+1).zfill(2)
            + '-ptl_' + petal_id + '-df_transformations.csv'))
    df_transformations.to_pickle(os.path.join(
            fig_save_dir,
            str(petal_ids.index(petal_id)+1).zfill(2)
            + '-ptl_' + petal_id + '-df_transformations.pickle'),
        compression='gzip')


def ics_output():

    # per DESI-2850
    df_focal_plane_ics = pd.DataFrame(
            index=np.arange(10),
            columns=['petal_id', 'petal_loc',
                     'petal_offset_x', 'petal_offset_y', 'petal_offset_z',
                     'petal_rot_x', 'petal_rot_y', 'petal_rot_z',
                     'tb0_543_x', 'tb0_543_y', 'tb0_543_z',
                     'tb1_544_x', 'tb1_544_y', 'tb1_544_z',
                     'tb2_545_x', 'tb2_545_y', 'tb2_545_z'])

    # convention is R = matmul(Rz(gamma), Ry(beta), Rx(alpha))
    for petal_location in petal_locations:

        # incorrect transformation values, needs work
        petal_id = petal_id_lookup[petal_location]
        df = pd.read_pickle(os.path.join(
                fig_save_dir,
                str(petal_ids.index(petal_id)+1).zfill(2)
                + '-ptl_' + petal_id + '-df_data.pickle'),
            compression='gzip')
        dtb_pos_ptl = np.array([df.loc[('dtb_x', 'ZBF'), 'actual'][0:3],
                                df.loc[('dtb_y', 'ZBF'), 'actual'][0:3],
                                df.loc[('dtb_z', 'ZBF'), 'actual'][0:3]])
        dtb_pos_act_cs5 = dtb_pos_cmm[petal_location].T

        # find transformation from petal local CS to focal plane CS5
        def total_residue_ics(parameters):
            # this is the function to be minimised
            alpha = parameters[0]
            beta = parameters[1]
            gamma = parameters[2]
            T = parameters[3:]
            R = Rxyz(alpha, beta, -gamma)  # yaw-pitch-roll system
            dtb_pos_rot = np.empty(dtb_pos_ptl.shape)
            for j in range(dtb_pos_ptl.shape[1]):
                dtb_pos_rot[:, j] = matmul(R, dtb_pos_ptl[:, j]) + T
            return np.sum(
                    np.linalg.norm(dtb_pos_rot - dtb_pos_act_cs5), axis=0)

        gamma0 = 2*np.pi/10*(petal_location-3)
        p0 = np.array([0.0, 0.0, gamma0, 0.0, 0.0, 0.0])
        bounds = ((-np.pi/2, np.pi/2), (-np.pi/2, np.pi/2),
                  (gamma0-1, gamma0+1),
                  (-10, 10), (-10, 10), (-10, 10))
        sol_ics = minimize(total_residue_ics, p0,
                           bounds=bounds, method='SLSQP',
                           options={'ftol': 1e-12,
                                    'eps': 1e-12,
                                    'disp': True,
                                    'maxiter': 10000})
        # reset angle within the range (-pi, pi)
        sol_ics.x[2] = np.arctan2(np.sin(sol_ics.x[2]),
                                  np.cos(sol_ics.x[2]))
        print('PTL{} CS5 transformation for ICS found as \n'.format(petal_id)
              + 'X Rotation (Roll) : {}° \n'.format(np.degrees(sol_ics.x[0]))
              + 'Y Rotation (Pitch): {}° \n'.format(np.degrees(sol_ics.x[1]))
              + 'Z Rotation (Yaw)  : {}° \n'.format(np.degrees(sol_ics.x[2]))
              + 'Translation: {} \n'.format(sol_ics.x[3:])
              + 'With least square: {} \n'.format(sol_ics.fun))

        df_focal_plane_ics.loc[petal_location] = {
                'petal_id':       petal_id_lookup[petal_location],
                'petal_loc':      petal_location,
                'petal_offset_x': sol_ics.x[3],
                'petal_offset_y': sol_ics.x[4],
                'petal_offset_z': sol_ics.x[5],
                'petal_rot_x':    np.degrees(sol_ics.x[0]),
                'petal_rot_y':    np.degrees(sol_ics.x[1]),
                'petal_rot_z':    np.degrees(sol_ics.x[2]),
                'tb0_543_x':      dtb_pos_cmm[petal_location][0, 0],
                'tb0_543_y':      dtb_pos_cmm[petal_location][0, 1],
                'tb0_543_z':      dtb_pos_cmm[petal_location][0, 2],
                'tb1_544_x':      dtb_pos_cmm[petal_location][1, 0],
                'tb1_544_y':      dtb_pos_cmm[petal_location][1, 1],
                'tb1_544_z':      dtb_pos_cmm[petal_location][1, 2],
                'tb2_545_x':      dtb_pos_cmm[petal_location][2, 0],
                'tb2_545_y':      dtb_pos_cmm[petal_location][2, 1],
                'tb2_545_z':      dtb_pos_cmm[petal_location][2, 2]}

    df_focal_plane_ics.to_csv(
            os.path.join(fig_save_dir, 'focal_plane_metrology_ics.csv'),
            columns=['petal_id', 'petal_loc',
                     'petal_offset_x', 'petal_offset_y', 'petal_offset_z',
                     'petal_rot_x', 'petal_rot_y', 'petal_rot_z',
                     'tb0_543_x', 'tb0_543_y', 'tb0_543_z',
                     'tb1_544_x', 'tb1_544_y', 'tb1_544_z',
                     'tb2_545_x', 'tb2_545_y', 'tb2_545_z'],
            float_format='%.15f')
    df_focal_plane_ics.to_pickle(
            os.path.join(fig_save_dir, 'focal_plane_metrology_ics.pickle'),
            compression='gzip')

# %% main


if __name__ == '__main__':

    # evaluate quality of petal-ring integration by calculating throughput
    pool = Pool()  # create a multiprocessing Pool
    pool.map(evaluate_petal, petal_locations)  # evaluate petal throughput
    ics_output()
