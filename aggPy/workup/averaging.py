#averaging.py
"""

    Average a provided data type from pandas df

"""

import numpy as np

#####

def avg_list(data):
    values = []
    for value in data:
        for v in value:
            values.append(v)
    return np.mean(values), np.std(values)

#####

def avg_dict(data):
    total = {}
    stdev = {}
    for value in data:
        for k,v in value.items():
            try:
                total[k].append(v)
            except KeyError:
                total[k] = [v]
    
    for k in total.keys():
        stdev[k] = np.std(total[k])
        total[k] = np.mean(total[k])

    return total, stdev

#####
