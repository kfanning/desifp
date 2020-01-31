import sys
sys.path.append("/data/focalplane/pecs_multi_petal/")
import numpy as np
import pandas as pd
from poscalibrationfits import PosCalibrationFits

df = pd.read_pickle('/n/home/desiobserver/enabled_pos.pkl.gz')
randposlocXY = np.random.rand(len(df)*10, 2) * 3.5
poslocR = np.linalg.norm(randposlocXY, axis=1)
mask = poslocR < 3.5
poslocXY = randposlocXY[mask]
df['poslocX'] = poslocXY[:len(df), 0]
df['poslocY'] = poslocXY[:len(df), 1]
fit = PosCalibrationFits(use_doslib=True)
fit.init_posmodels(posids=df.index)
for posid, row in df.iterrows():
	trans = fit.posmodels[posid].trans
	QS = trans.poslocXY_to_QS(row[['poslocX', 'poslocY']])
	obsXY = trans.QS_to_obsXY(QS, cast=True)
	df.loc[posid, 'Q'], df.loc[posid, 'S'] = QS[0], QS[1]
	df.loc[posid, 'obsX'], df.loc[posid, 'obsY'] = obsXY[0], obsXY[1]
df.to_pickle('rand_targets.pkl.gz')
df.to_csv('rand_targets.csv')
