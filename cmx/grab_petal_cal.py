import psycopg2
from astropy.table import Table
from astropy.stats import sigma_clipped_stats
import datetime
from pandas.io.sql import read_sql

def get_petal_calib(petal_loc, firstlast='last',pandas=True, cols=''):
    '''
    returns Table of calibration data for the requested petal location
    pandas = True returns pandas df. else astropy table
    '''
    #- which PETAL_ID is plugged into which PETAL_LOC?
    #- petal_loc2id[loc] = id
    #- petal_id2loc[id] = loc
    #- From DESI-5286
    petal_loc2id = {2:6, 9:9, 0:4, 3:3, 8:7, 4:8, 6:11, 7:2, 5:10, 1:5}
    petal_id2loc = dict()
    for ploc, pid in petal_loc2id.items():
        petal_id2loc[pid] = ploc
    
    petal_id = petal_loc2id[petal_loc]
    print('Getting latest calib for petal_id={} installed at petal_loc={}'.format(petal_id, petal_loc))
    
    if firstlast == 'last':
        minmax = 'max'
    elif firstlast == 'first':
        minmax = 'min'
    else:
        raise ValueError('firstlast should be "first" or "last"')
    if not(cols):
        cols = 't.length_r1, t.length_r2, t.offset_x, t.offset_y, t.offset_t, t.offset_p'
    query = """
    select t.petal_id, t.device_loc, t.pos_id, t.time_recorded,
        {cols}
    from posmovedb.positioner_calibration_p{petal_id} t
    inner join (
        select device_loc, {minmax}(time_recorded) as MaxDate
        from posmovedb.positioner_calibration_p{petal_id}
        group by device_loc
    ) tm on t.device_loc = tm.device_loc and t.time_recorded = tm.MaxDate
    """.format(cols=cols, petal_id = petal_id, minmax=minmax)

    comm = psycopg2.connect(host='db.replicator.dev-cattle.stable.spin.nersc.org', port=60042,
                           database='desi_dev', user='desi_reader', password='reader')

    cx = comm.cursor()
    cx.execute(query)
    results = cx.fetchall()

    if pandas:
        tx = read_sql(query, comm)
    else:
        tx = Table(rows=results,
                   names=['petal_id', 'device_loc', 'pos_id', 'time_recorded',
                          'r1', 'r2', 'x', 'y'],
                   dtype=[int, int, str, datetime.datetime, float, float, float, float])

        tx.meta['PETAL_ID'] = petal_id
        tx.meta['PETAL_LOC'] = petal_loc

        tx.sort('device_loc')

        tmin = np.min(tx['time_recorded'])
        tmax = np.max(tx['time_recorded'])
        print('includes calibrations uploaded from {} to {}'.format(
            tmin.strftime('%Y-%m-%d %H:%M:%S'),
            tmax.strftime('%Y-%m-%d %H:%M:%S'),
        ))
    return tx

# Try it out
ptl0_calib = get_petal_calib(0)