import numpy as np
import json 
import random
import os

from pathlib import Path

####
def read_file(filename):
    current_dir = Path(__file__).parent

    file_path = current_dir / filename

    with open(file_path, 'r', encoding='utf-8') as f:
        file = [line.strip().split() for line in f]
        return file, file_path

#####
def write_atomtypes(lt_file, p=''):
    types = []
    print('#####')
    print(p, lt_file)
    try:
        with open(f'{p}/{lt_file}', 'r') as f:
            lines = [line.strip().split() for line in f]
    except FileNotFoundError:
        try:
            lines, file_path = read_file(f'mol_templates/{lt_file}')
            os.system(f'cp {file_path} {p}')
        except FileNotFoundError:
            print(f'{lt_file} not found')

    for line in lines:
        for l in line:
            try:
                if l != '#' and l.split(':')[0] == '@atom':
                    types.append(l.split(':')[1])
            except IndexError:
                pass

    return list(set(types))

#####
def make_metaAtomTypes(solvents, molecules, path, output='meta_atomTypes.json'):
    types = {}
    with open(f'{path}/{output}','w') as f:
        for s, solvent in enumerate(solvents):
            types[f'solvent_{s}'] = solvent.split(".")[0]
            atom_types = write_atomtypes(solvent, p=path)
            #[f.write(f'{atom} - solvent_{s}\n') for atom in atom_types]
            for atom in atom_types:
                try: types[f'{solvent.split(".")[0]}'].append(atom)
                except KeyError: types[f'{solvent.split(".")[0]}'] = [atom]
        for m, molecule in enumerate(molecules):
            types[f'molecule_{m}'] = molecule.split(".")[0]
            atom_types = write_atomtypes(molecule, p=path)
            #[f.write(f'{atom} - molecule_{m}\n') for atom in atom_types]
            for atom in atom_types:
                try: types[f'{molecule.split(".")[0]}'].append(atom)
                except KeyError: types[f'{molecule.split(".")[0]}'] = [atom] 


        json.dump(types, f, indent=4)

#####
def make_grid(n_mols, dim):
    x = np.linspace(0, dim, n_mols+1, endpoint=False)[1:]

    x = np.array(list(zip(x,x,x)))

    grid = np.zeros((n_mols, 3))
    for i in range(0,3):
        coord = np.random.choice(x[:,i], size=n_mols, replace=False)
        grid[:,i] = coord

    return grid

#####

def make_box(path, molecules=[], solvents=[], n_molecules=[0], n_solvents=[0], box_side=50):

    make_metaAtomTypes(solvents, molecules, path)

    with open(f"{path}/system.lt","w") as f:
        f.write('# -- System --\n\n')
        for solvent in solvents:
            f.write(f'import "./{solvent}"\n')
        for molecule in molecules:
            f.write(f'import "./{molecule}"\n')
        f.write('\n')
        
        f.write('write_once("Data Boundary") {\n')
        f.write(f'\t 0 {box_side} xlo xhi\n')
        f.write(f'\t 0 {box_side} ylo yhi\n')
        f.write(f'\t 0 {box_side} zlo zhi\n')
        f.write('}\n\n')

        molecules = solvents + molecules
        n_molecules = n_solvents + n_molecules
        print(n_molecules)
        for i, n_mol in enumerate(n_molecules):
            if n_mol == 0: continue
            name = molecules[i].split('.')[0]
            grid = make_grid(n_mol, box_side)
            for _ in range(0, n_mol):
                coords = grid[_]
                theta = round(random.uniform(0, 360))
                axis = [random.randint(0, 1) for a in range(0, 3)]
                f.write(f'{name}{_} = new {name}')
                f.write(f'.rot({theta}, {axis[0]}, {axis[1]}, {axis[2]})')
                f.write(f'.move({coords[0]}, {coords[1]}, {coords[2]})\n')


#####
if __name__ == '__main__':
    make_box(
            molecules=['myMol.lt'],
            solvents=['spce.lt'],
            n_molecules=[100],
            n_solvents=[0],
            box_side=50,
        )

