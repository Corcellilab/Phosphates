#dump_analyze_csv.py

import pandas as pd

from .averaging import *

####
def dump_analyze_csvs(results, outfile):
    results.to_csv(f'{outfile}.csv', index=True)

    avg_results = {}
    for i in results.columns:
        typing = type(results[i][0])
        if typing == list:
            avg, stdev = avg_list(results[i][:])
        elif typing == dict:
            avg, stdev = avg_dict(results[i][:])
        else:
            numeric_values = [x for x in results[i] if isinstance(x, (int, float))]
            avg = pd.Series(numeric_values).mean()
            stdev = pd.Series(numeric_values).std()

        avg_results[f'avg_{i}'] = avg
        avg_results[f'stdev_{i}'] = stdev

    pd.DataFrame(avg_results).sort_index().to_csv(f'avg_{outfile}.csv')

####

