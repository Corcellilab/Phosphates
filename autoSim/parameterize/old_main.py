from .make_molecule import make_molecule
from .make_gcrt import make_gcrt
from .edit_submit import edit_submit
from .make_box import make_box
from .edit_frcmod import edit_frcmod

import sys
import os
import json
import subprocess
from pathlib import Path

####
def read_file(filename):
    current_dir = Path(__file__).parent

    file_path = current_dir / filename

    return file_path

######
def main_Parm():

    root = os.getcwd()

    p = Path(f'{root}/parm')
    try:
        p.mkdir(exist_ok=False)
    except FileExistsError:
        os.system(f'rm -r parm')
        p.mkdir()

    with open(f"{root}/settings.json",'r') as j:
        settings = json.load(j)

    for molecule in settings['parm']['smiles']:
    smiles = settings['parm']['smiles']
    name = settings['parm']['name']

    print('Making xyz...')
    xyz_file = make_molecule(smiles, output=f'{p}/{name}')
    print('Done')

    #Get charge returned from  make_molecules?#
    master_gcrt = 'master_gcrt'
    gcrt_file = make_gcrt(
            xyz_file=xyz_file, 
            master_gcrt=master_gcrt, 
            output=f'{p}/{name}',
            charge=settings['parm']['charge'],
            spin=settings['parm']['spin'],
        )

    print('Run parameterization...')
    print('\t run gauss...')
    os.system(f'g16 {p}/{name}.gcrt')
    print('\t done')
    print('\t run antechamber...')
    os.system(f'antechamber -i {p}/{name}.gesp -fi gesp -o {p}/{name}.mol2 -fo mol2 -c resp -pl 30 -s 2')
    print('\t done')
    print('\t parmchk...')
    os.system(f'parmchk2 -i {p}/{name}.mol2 -f mol2 -o {p}/{name}.frcmod')
    print('\t done')
    print('Done')

    #Build system
    edit_frcmod(frcmod=f'{p}/{name}.frcmod')

    ff_path = read_file('mol_templates/gaff2.lt')

    #os.system(f'python3 amber2lt.py --in {p}/{name}.frcmod --name gaff2')
    #os.system(f'python3 mol22lt.py --in {p}/{name}.mol2 --out {p}/{name}.lt --name {name} --ff gaff2 --ff-file gaff2.lt')

    from .amber2lt import main as amber2lt_main
    amber2lt_main(f'{p}/{name}.frcmod')

    from .mol22lt import main as mol22lt_main
    mol22lt_main(f'{p}/{name}.mol2', f'{p}/{name}.lt', f'{name}', 'gaff2', f'{ff_path}')

    os.system(f'mv A* {name}.* parm')
    os.system(f'mv QOUT parm')
    os.system(f'mv *out* parm')
    os.system(f'mv punch parm')

    os.system(f'ln -sv {ff_path} {p}')

    #if settings['parm']['charge'] > 0:
    #   molecules = [f'{name}.lt', 'cl.lt']
    #elif settings['parm']['charge'] < 0:
    #   molecules = [f'{name}.lt', 'na.lt']

    print('Build system...')
    make_box(path=p,
            molecules=[f'{name}.lt'],
            solvents=['spce.lt'],
            n_molecules=settings['parm']['n_molecules'],
            n_solvents=settings['parm']['n_solvents'],
            box_side=settings['parm']['box_side']
        )

    os.system(f'moltemplate.sh -atomstyle full {p}/system.lt')

    os.system(f"rm {root}/run.in.EXAMPLE")

    print('Done')

    print('Done parameterization')

#####
if __name__ == "__main__":
    #smiles = "CCl"
    #main_Parm(smiles=smiles, name='c')
    main_Parm()

