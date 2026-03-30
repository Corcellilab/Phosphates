from read_lmp_dump import read_lmp_dump
import json
import os

#####
def run_defrost(run=400000, simTemp=300, simPress=300):
    with open('in.defrost','r') as f:
        lines = [line.strip().split() for line in f]
    
    for line in lines:
        try:
            if line[0] == 'variable':
                if line[1] == 'run':
                    line[-1] == f'{run}'
                elif line[1] == 'simTemp':
                    line[-1] == f'{simTemp}'
                elif line[1] == 'simPress':
                    line[-1] == f'{simPress}'
        except IndexError:
            pass

    with open('in.defrost','w') as f:
        for line in lines:
            f.write(' '.join(line))
            f.write('\n')

    os.system('lmp -in in.defrost')

    filename = 'defrost.txt'
    prefix = filename.split(".")[0]

    data = read_lmp_dump(filename)
    avgs = graph_lmp_dump(data, window_size=10, outfile=f'{prefix}_avg.png')

    with open(f'{prefix}.json','w') as j:
        json.dump(avgs, j, indent=4) 

    return avgs

#####

if __name__ == '__main__':
    avgs = run_defrost(run=100, simTemp=300, simPress=300)
