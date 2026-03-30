from .run_prod import run_prod
from .find_port import find_port
from .velOpt import main_vel

import json
import os
from pathlib import Path

#####

def read_file(filename):
    current_dir = Path(__file__).parent

    file_path = current_dir / filename

    return file_path

######

def main_Prod(settings='settings.json'):

    root = os.getcwd()

    with open(f'{root}/{settings}','r') as f:
        settings = json.load(f)

    p = Path(f'{root}/prod')
    try:
        p.mkdir()
    except FileExistsError:
        os.system(f'rm -r prod')
        p.mkdir()

    lmp_ins_p = read_file('lmp_ins') 
    
    os.system(f'ln -sv {lmp_ins_p} {p}')

    lmp_p = f'{p}/lmp_ins'

    os.system(f'cp {root}/equib/nve/final.preNVE.data {root}/prod')
    os.system(f'cp {lmp_p}/in.prod {p}')

    os.chdir('prod')

    #Get solute data file
    if settings['GPU'] == 'False': os.system(f'lmp -in {lmp_p}/in.delete')
    elif settings['GPU'] == 'True': os.system(f'lmp -in {lmp_p}/in.delete -sf gpu')

    #Update with Port number#
    port = find_port()
  
    print(port)

    run_prod(lmp_file=f'{p}/in.prod',
        run=settings['prod']['run'],
        simTemp=settings['prod']['temp'],
        simPress=settings['prod']['press'],
        port=port,
    )

    #Run prod in background
    #Connect to simulation client with main_vel()
    os.system(f'hostname -i > ip.txt')

    if settings['GPU'] == 'False': os.system(f"mpirun -np 1 lmp_mpi -in {p}/in.prod &")
    elif settings['GPU'] == 'True': os.system(f"mpirun -np 1 lmp_mpi -in {p}/in.prod -sf gpu &")

    os.system(f'python3 {lmp_p}/velOpt.py -r {root}')

    #main_vel(
    #        data='solute.data',
    #        max_dt=settings['spectra']['max_dt'],
    #        max_t=settings['spectra']['max_t'],
    #        jump_t=1,
    #        dt=0.5e-15
    #    )

#####

if __name__ == '__main__':
    settings = input('Settings JSON file: ')
    main_Prod(settings)

