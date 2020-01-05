''' this finds the pos ids of extra spots in fvc image'''
import os
import numpy as np
import pandas as pd
from petaltransforms import PetalTransforms
from DOSlib.constants import ConstantsDB
import matplotlib.pyplot as plt

# get pos file xy columns saved for ds9 overlay
# but with the y flip and petal rotation is hard to inspect
pcid = 6
path = "/data/focalplane/logs/kpno/20200104/00037869/fvc_proxy_centers_20200105011630060034.json"
c = pd.read_json(path).T
c[['x', 'y']] += 3000
pos = c[c['comment'] == 'fiber']
pos.to_csv(os.path.join(os.path.dirname(path), f'pc{pcid:02}-matched.pos'), index=False, sep='\t')

# compare nonintact in DB with positioner index table
path = os.getenv('DOS_POSITIONERINDEXTABLE')
if not os.path.isfile(path):
    raise ValueError(f'$DOS_POSITIONERINDEXTABLE = {path}')
pi_df = pd.read_csv(path)
pc = pi_df[(pi_df['PETAL_LOC']==7) & (pi_df['DEVICE_TYPE']=='POS')].set_index('DEVICE_ID')
# get all fibre non intact positioners from constantsDB
group, tag, snapshot = 'fiber_positioner', 'CURRENT', 'DESI'
df = pd.DataFrame.from_dict(ConstantsDB().get_constants( 
	group=group, tag=tag, snapshot=snapshot)[group], orient='index')
nonintact = df[df['FIBER_INTACT']==False]  # all nonintact
nonintact = nonintact[nonintact.index.isin(pc.index)]  # filter out specified pc
nonintact = pc.loc[nonintact.index]
print(nonintact)

# plot nonintact positioners
path = os.path.join(os.getenv('PLATE_CONTROL_DIR',
                              '/software/products/plate_control-trunk'),
                    'petal', 'positioner_locations_0530v14.csv')
ptlXYZ = pd.read_csv(path, usecols=['device_location_id', 'X', 'Y', 'Z'],
		             index_col='device_location_id')
ptlXYZ.index.rename('DEVICE_LOC', inplace=True)
trans = PetalTransforms(gamma=np.pi/5*(pcid-3))
obsXY = trans.ptlXYZ_to_obsXYZ(ptlXYZ.values.T).T[:, :2]
xy = pd.DataFrame(obsXY, index=ptlXYZ.index, columns=['obsX', 'obsY'])
xy['fvcX'] = -xy['obsX']
xy['fvcY'] = xy['obsY']
nonintact = nonintact.merge(xy, how='left', left_on='DEVICE_LOC', right_index=True) 
fig, ax = plt.subplots()
ax.plot(xy['fvcX'], xy['fvcY'], 'o')
ax.plot(nonintact['fvcX'], nonintact['fvcY'], 'rx')
for device_id, row in nonintact.iterrows():
	ax.annotate(f'{device_id} ({row.DEVICE_LOC})', xy=(row.fvcX, row.fvcY), 
		        xytext=(0, -5), textcoords='offset points')
ax.set_aspect('equal')
fig.savefig(f'pc{pcid:02}_nonintact.pdf')

