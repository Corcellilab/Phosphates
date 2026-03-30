import os
import sys


def rewrite_settings(root):

    os.system(f'cp {root}/system.in.settings {root}/equib/backup_settings')

    with open(f'{root}/system.in.settings','r') as f:
        lines = [line.strip().split() for line in f]

    with open(f'{root}/system.in.settings','w') as f:
        for i,line in enumerate(lines):
            if 'harmonic' in line:
                del line[line.index('harmonic')]
            elif 'lj/charmm/coul/long' in line:
                del line[line.index('lj/charmm/coul/long')]
            elif 'fourier' in line:
                del line[line.index('fourier')]
            elif 'cvff' in line:
                del line[line.index('cvff')]

            f.write(" ".join(map(str, line)) + '\n')

######
if __name__ == '__main__':
    rewrite_settings(os.getcwd())
