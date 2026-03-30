from .read_lmp_dump import read_lmp_dump
from .graph_lmp_dump import graph_lmp_dump

from random import randint

import json
import os

#####
def run_lmp(lmp_file='', gpu=False, run=400000, simTemp=300, simPress=1, window_size=10, box_side=50):
    with open(lmp_file,'r') as f:
        lines = [line.strip().split() for line in f]
    
    for line in lines:
        try:
            if line[0] == 'variable':
                if line[1] == 'run':
                    line[-1] = f'{run}'
                elif line[1] == 'simTemp':
                    line[-1] = f'{simTemp}'
                elif line[1] == 'simPress':
                    line[-1] = f'{simPress}'
                elif line[1] == 'rnseed':
                    line[-1] = f'{randint(0, 100)}'
                elif line[1] == 'boxSide':
                    line[-1] = f'{box_side}'
        except IndexError:
            pass

    with open(lmp_file,'w') as f:
        for line in lines:
            f.write(' '.join(line))
            f.write('\n')
    
    print('################################')
    print(f'Running: {lmp_file}')
    if gpu == 'False': os.system(f'lmp -in {lmp_file}')
    elif gpu == 'True': os.system(f'lmp -in {lmp_file} -sf gpu')
    print('Done')

    filename = f"{lmp_file.split('.')[1]}.txt"
    prefix = filename.split(".")[0]

    data = read_lmp_dump(filename)
    avgs = graph_lmp_dump(data, window_size=window_size, outfile=f'{prefix}_avg.png')

    with open(f'{prefix}.json','w') as j:
        json.dump(avgs, j, indent=4) 

    os.system('cp final.*.data ../')

    return avgs

#####

if __name__ == '__main__':
    avgs = run_defrost(run=100, simTemp=300, simPress=300)
