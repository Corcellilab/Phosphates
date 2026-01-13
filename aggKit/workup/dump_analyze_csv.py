#dump_analyze_csv.py

import pandas as pd

from .averaging import *


####
def dump_analyze_csvs(r, outfile):
    print('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')
    r.to_csv(f'{outfile}.csv', index=True)

    avg_results = {}
    for i in r.columns:
        typing = type(r[i][0])
        if typing == list:
            avg, stdev = avg_list(r[i][:])
        elif typing == dict:
            avg, stdev = avg_dict(r[i][:])
        else:
            numeric_values = [x for x in r[i] if isinstance(x, (int, float))]
            avg = pd.Series(numeric_values).mean()
            stdev = pd.Series(numeric_values).std()

        avg_results[f'avg_{i}'] = avg
        avg_results[f'stdev_{i}'] = stdev


    result = pd.DataFrame(avg_results)
   
    result.sort_index().to_csv(f'avg_{outfile}.csv')

    return avg_results

####
