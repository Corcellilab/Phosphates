from .read_lmp_dump import read_lmp_dump
from .graph_lmp_dump import graph_lmp_dump
from .find_port import find_port

from random import randint

import json
import os

#####
def run_prod(lmp_file='', run=400000, simTemp=300, simPress=1, window_size=10, port=0,):
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
                elif line[1] == 'port':
                    line[-1] = f'{port}'
        except IndexError:
            pass

    with open(lmp_file,'w') as f:
        for line in lines:
            f.write(' '.join(line))
            f.write('\n')

#####

if __name__ == '__main__':
    
    settings = 'settings.json'

    with open(settings,'r') as f:
        settings = json.load(f)
    
    port = find_port()

    run_prod(lmp_file='in.prod',
        run=settings['prod']['run'],
        simTemp=settings['prod']['temp'],
        simPress=settings['prod']['press'],
        port = port,
    )

