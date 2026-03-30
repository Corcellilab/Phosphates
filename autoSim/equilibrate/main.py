from .run_lmp import run_lmp
from .check_convg import check_convg
from .edit_min import edit_min
from .rewrite_init import rewrite_init
from .rewrite_settings import rewrite_settings

import os
import json
from pathlib import Path

#####

def read_file(filename):
    current_dir = Path(__file__).parent

    file_path = current_dir / filename

    return file_path

#####
def not_convg_error(col, avg, stdev):
    print(f'WARNING: {col} not converged, average value: {avg} +/- {2*stdev}')

#####
def main_Equib(t_temp=300, t_press=1, settings='settings.json'):

    root = os.getcwd()

    with open(f'{root}/{settings}','r') as f:
        settings = json.load(f)

    if settings['GPU'] == 'True':
        rewrite_init(root)
        rewrite_settings(root)

    p = Path(f'{root}/equib')
    try:
        p.mkdir(exist_ok=False)
    except FileExistsError:
        os.system(f'rm -r equib')
        p.mkdir()
   
    lmp_ins_p = read_file('lmp_ins')

    os.system(f'ln -sv {lmp_ins_p} {p}')

    lmp_p = f'{p}/lmp_ins'

    os.chdir('equib')

    #Minimize#
    edit_min(root)
    print('Run minimization...')
    os.system('mkdir minimization')
    os.chdir('minimization')
    print(settings)
    if settings['GPU'] == 'False': os.system(f'lmp -in {lmp_p}/in.min')
    elif settings['GPU'] == 'True': os.system(f'lmp -in {lmp_p}/in.min -sf gpu')
    print('Done')
    os.chdir('../')

    
    #Defrost#
    os.system('mkdir defrost')
    os.chdir('defrost')
    defrost_avgs = run_lmp(lmp_file=f'{lmp_p}/in.defrost', 
            run=settings['defrost']['run'], 
            simTemp=settings['defrost']['temp'], 
            simPress=settings['defrost']['press'],
            gpu=settings['GPU']
        )
    os.chdir('../')

    #NPT#
    os.system('mkdir npt')
    os.chdir('npt')
    npt_avgs = run_lmp(lmp_file=f'{lmp_p}/in.eqNPT', 
            run=settings['npt']['run'], 
            simTemp=settings['npt']['temp'], 
            simPress=settings['npt']['press'],
            gpu=settings['GPU']
        )

    convg, avg, stdev = check_convg('eqNPT.json', 'Temp', start=1/3, t_value=t_temp)
    if convg == False: not_convg_error('Temp', avg, stdev)
    if avg <= 10: print('WARNING: Box length shorter than cutoff distance')
    
    #check output txt and pipe box side into eqNVT
    convg, avg, stdev = check_convg('eqNPT.json', 'Volume', start=1/2, t_value='avg')
    if convg == False: not_convg_error('Volume', avg, stdev)

    os.chdir('../')

    #NVT#
    os.system('mkdir nvt')
    os.chdir('nvt')
    nvt_avgs = run_lmp(lmp_file=f'{lmp_p}/in.eqNVT', 
            run=settings['nvt']['run'], 
            simTemp=settings['nvt']['temp'], 
            simPress=settings['nvt']['press'], 
            box_side=avg**(1/3),
            gpu=settings['GPU']
        )
    convg, avg, stdev = check_convg('eqNVT.json', 'Temp', start=1/3, t_value=t_temp)
    if convg == False: not_convg_error('Temp', avg, stdev)

    os.chdir('../')

    #NVE#
    os.system('mkdir nve')
    os.chdir('nve')
    nve_avgs = run_lmp(lmp_file=f'{lmp_p}/in.preNVE', 
            run=settings['nve']['run'], 
            simTemp=settings['nve']['temp'], 
            simPress=settings['nve']['press'],
            gpu=settings['GPU']
        )
    convg, avg, stdev = check_convg('preNVE.json', 'Temp', start=1/3, t_value=t_temp)
    if convg == False: not_convg_error('Temp', avg, stdev)

    os.chdir('../')
    os.chdir(f'{root}')

#####
if __name__ == '__main__':

    #Run args
    #run_time - simTemp - simPress 

    t_temp = 1
    t_press = 1

    defrost_avgs = run_lmp(lmp_file='in.defrost', run=1000, simTemp=t_temp, simPress=t_press)

    npt_avgs = run_lmp(lmp_file='in.eqNPT', run=1000, simTemp=t_temp, simPress=t_press)

    convg, avg, stdev = check_convg('eqNPT.json', 'Temp', start=1/3, t_value=t_temp)
    if convg == False: not_convg_error('Temp', avg, stdev)

    #check output txt and pipe box side into eqNVT
    convg, avg, stdev = check_convg('eqNPT.json', 'Volume', start=1/2, t_value='avg')
    if convg == False: not_convg_error('Volume', avg, stdev)

    nvt_avgs = run_lmp(lmp_file='in.eqNVT', run=1000, simTemp=t_temp, simPress=t_press, box_side=avg**(1/3))

    nve_avgs = run_lmp(lmp_file='in.preNVE', run=1000, simTemp=t_temp, simPress=t_press,)

