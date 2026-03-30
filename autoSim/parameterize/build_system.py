from edit_frcmod import edit_frcmod
from make_box import make_box

import os

def main(name='c'):

    edit_frcmod(frcmod=f'{name}.frcmod')

    os.system(f'python3 amber2lt.py --in {name}.frcmod --name gaff2')
    os.system(f'python3 mol22lt.py --in {name}.mol2 --out {name}.lt --name {name} --ff gaff2 --ff-file gaff2.lt')

    make_box(
            molecules=[f'{name}.lt'],
            solvents=['spce.lt'],
            n_molecules=[100],
            n_solvents=[1000],
            box_side=40
        )

    os.system('moltemplate.sh -atomstyle full system.lt')

#####

if __name__ == '__main__':
    main()

