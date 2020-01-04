import pandas as pd
from DOSlib.proxies import Petal

path = '/data/focalplane/logs/calib_logs/20191220T150238-0700-arc_calibration.pkl'
df = pd.read_pickle(path)
for i in range(10):
    ptl = Petal(i)
    ptl.set_calibration(df, tag='OLD')

