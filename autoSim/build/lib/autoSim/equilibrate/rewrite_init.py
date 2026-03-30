import sys
import os

def rewrite_init(root):
    
    os.system(f'cp {root}/system.in.init {root}/equib/backup_system.in.init')

    with open(f'{root}/system.in.init','r') as f:
        lines = [line.strip().split() for line in f]

    with open(f'{root}/system.in.init','w') as f:
        for i,line in enumerate(lines):
            try:
                ind = line.index('hybrid')
                del line[ind]
            except ValueError:
                pass

            f.write(" ".join(map(str, line)) + '\n')

#####
if __name__ == '__main__':
    root = os.getcwd()
    rewrite_init(root)

