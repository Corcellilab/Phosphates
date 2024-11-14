#time_corr.py
"""

    Calculate autocorrelation function for
    dictonary of {atomic index: [(x1,y1,z1),(x2,y2,z2),...]}

"""

import numpy as np
import matplotlib.pyplot as plt
import sys

def time_corr_vector(data, min_dt=0, max_dt=50, skip_dt=1, t_prop=1):

    max_dt = int(max_dt/5)

    keys = [_ for _ in range(min_dt, max_dt, skip_dt)]
    values = [[] for _ in range(len(keys))]
    values = dict(zip(keys, values))

    velocities = [vector / np.linalg.norm(vector) for vector in data]

    print('Time Correlation Calculating...')

    for i in range(min_dt, max_dt, t_prop):
        
        v1 = velocities[i]
        for j in range(min_dt, max_dt, skip_dt):
            try: v2 = velocities[i+j]
            except IndexError: break
            values[j].append(np.dot(v1,v2))

    x = list()
    y = list()
    for k,v in values.items():
        x.append(np.average(k))
        y.append(np.average(v))

    return x, y


